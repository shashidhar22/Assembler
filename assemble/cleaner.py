import os
import re
import csv
import sys
import Bio
import glob
import time
import timeit
import logging
import argparse
import itertools
import subprocess
import collections
import configparser
from Bio import SeqIO
from itertools import repeat
from multiprocessing import Pool
from collections import defaultdict

#All pacbio features will be reintroduced in version 0.9.1
logger = logging.getLogger('Nutcracker')

class Cleaner:
    def __init__(self, bowtie_path, bbduk_path, out_path, threads):
        
        logger = logging.getLogger('Nutcracker')
        #Initialize values 
        self.bowtie_path = os.path.abspath(bowtie_path)
        self.bbduk_path = os.path.abspath(bbduk_path)
        self.out_path = os.path.abspath(out_path)
        self.threads = threads
        
        #Create output directory
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)
        return
        
    def buildIndex(self, reference):
        '''Build bowtie index if it doesnt exist'''
        
        #Initialize values
        reference = os.path.abspath(reference)
        btindex = os.path.splitext(reference)[0]


        #Check for index
        logger.info('Checking for index')
        conindex = glob.glob('{0}*bt2'.format(btindex))
        if conindex:
            logger.info('Bowtie index exists : {0}; Skipping step;'.format(btindex))
            return(0, btindex)
        
        #Prepare run command
        logger.info('Creating bowtie index')
        #Setup bowtie index command
        bicmd = ['{0}-build'.format(self.bowtie_path ), reference, btindex]
        logger.debug('Running bowtie index with the following command')
        logger.debug('{0}'.format(' '.join(bicmd)))
    
        #Running bowtie build
        birun = subprocess.Popen(bicmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        #Capture stdout and stderr
        bilog = birun.communicate()
        elapsed = timeit.default_timer() - start
        for logs in bilog:
            log = logs.decode('utf-8').split('\n')
            for lines in log:
                logger.debug(lines)
        if birun.returncode != 0 :
            logger.error('Bowtie index failed with exit code : {0}; Check runtime log for details.'.format(birun.returncode))
            return(birun.returncode, None)
        else:
            logger.info('Bowtie index completed successfully; Runtime : {0}'.format(elapsed))
            logger.info('Reference index can be found at : {0}.*.bt2'.format(btindex))
            logger.close()
            return(birun.returncode, btindex)
    
        
    def deconIllumina(self, read_one, read_two, btindex):
        '''Remove contaminants from illumina reads'''
        #Start logging
        start = timeit.default_timer() 
        
        #Prepare run commands
        cread_base = '{0}/sample_cleaned.fq'.format(self.out_path)
        cread_one = '{0}/sample_cleaned.1.fq'.format(self.out_path)
        cread_two = '{0}/sample_cleaned.2.fq'.format(self.out_path)
        cread_sam = '{0}/sample_aligned.sam'.format(self.out_path)
        logger.info('Decontamination  of illumina reads  started')
        
        #Setup bowtie command
        if read_two:
            bcmd = [self.bowtie_path, '-p', self.threads, '-x', btindex, '-1', read_one, '-2', read_two,
                    '--un-conc', cread_base, '-S', cread_sam]
            logger.debug('Running Bowtie with the following command')
            logger.debug('{0}\n'.format(' '.join(bcmd)))
        else:
            bcmd = [self.bowtie_path, '-p', self.threads, '-x', btindex, '-U', ','.join(read_one)
                    '--un-conc', cread_base, '-S', cread_sam]
            logger.debug('Running Bowtie with the following command')
            logger.debug('{0}\n'.format(' '.join(bcmd)))


    
        #Running bowtie
        brun = subprocess.Popen(bcmd, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, shell=False)
        
        #Capture stdout and stderr
        blog = brun.communicate()
        elapsed = timeit.default_timer() - start
        for logs in blog:
            log = logs.decode('utf-8').split('\n')
            for lines in log:
                logger.debug(lines)

    
        if brun.returncode != 0 :
            logger.error('Bowtie failed with exit code : {0}; Check runtime log for details.'.format(brun.returncode))
            return(brun.returncode, None, None, None)
        elif read_two:
            logger.info('Bowtie completed successfully; Runtime : {0}'.format(elapsed))
            logger.info('Decomtaminated reads can be found in : {0}'.format(self.out_path))
            return(brun.returncode, cread_one, cread_two, cread_sam)
        else:
            logger.info('Bowtie completed successfully; Runtime : {0}'.format(elapsed))
            logger.info('Decontaminated reads can be found in : {0}'.format(self.out_path))
            return(brun.returncode, cread_base, None, cread_sam)

    def trimIllumina(self, cread_one, cread_two, adapters):
        '''Trim adapters from illumina reads'''
        #Start logging
        start = timeit.default_timer()

        #Initialize values
        
        if cread_two:
            #If paired end
            tread_one = '{0}/cleaned_trimmed_r1.fastq'.format(self.out_path)
            tread_two = '{0}/cleaned_trimmed_r2.fastq'.format(self.out_path)
        else:
            #If single reads
            tread_one = '{0}/cleaned_trimmed.fastq'.format(self.out_path)
            tread_two = None

        #Prepare run commands
        logger.info('Trimming illumina reads\n')
        #Setup bbduk command
        if cread_two:
            tcmd = [self.bbduk_path, 'ref={0}'.format(','.join(adapters)), 'in={0}'.format(cread_one), 
                    'in2={0}'.format(cread_two), 'out1={0}'.format(tread_one), 
                    'out2={0}'.format(tread_two), 'ktrim=r', 'ktrim=1', 'k=27', 'mink=11', 'qtrim=rl',
                    'trimq=30', 'minlength=80', 'overwrite=t']
        else:
            tcmd = [self.bbduk_path, 'ref={0}'.format(','.join(adapters)), 'in={0}'.format(cread_one), 
                    'out1={0}'.format(tread_one), 'ktrim=r', 'ktrim=1', 'k=27', 'mink=11', 'qtrim=rl',
                    'trimq=30', 'minlength=80', 'overwrite=t']                                

        logger.debug('Running BBDuk with the following command:\n')
        logger.debug('{0}\n'.format(' '.join(tcmd)))

        #Runnig BBDuk
        trun = subprocess.Popen(tcmd, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, shell=False)

        #Capture the stdout and stderr
        tlog = trun.communicate()
        elapsed = timeit.default_timer() - start
        
    
        if trun.returncode != 0 :
            logger.error('BBDuk failed with exit code : {0}; Check runtime log for details.\n'.format(trun.returncode))
            return(trun.returncode, None, None)
        elif cread_two:
            logger.info('BBDuk completed successfully; Runtime : {0}\n'.format(elapsed))
            logger.info('Decontaminated reads can be found in : {0}\n'.format(self.out_path))
            return(trun.returncode, tread_one, tread_two)            
        else:
            logger.info('BBDuk completed successfully; Runtime : {0}\n'.format(elapsed))
            logger.info('Decontaminated reads can be found in : {0}\n'.format(self.out_path))
            return(trun.returncode, tread_one, None)


    def cleanSam(self, cread_sam):
        '''Remove sam file'''
        os.remove(cread_sam)
        return



if __name__ == '__main__':
    logger = logging.getLogger('Assembler')
    logger.setLevel(logging.DEBUG)

    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)

    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # add formatter to ch
    ch.setFormatter(formatter)

    # add ch to logger
    logger.addHandler(ch)

    logger.info('Starting assembly pipeline')
    logger.info('Running the following modes of analysis : Cleaner')

    bowtie_path = sys.argv[1]
    bbduk_path = sys.argv[2]
    prep_out_path = sys.argv[3]
    threads = '4'
    read_one = sys.argv[4]
    read_two = sys.argv[5]
    adapters = sys.argv[6].split(',')
    ref_path = sys.argv[7]
    cleaner = Cleaner(bowtie_path, bbduk_path, prep_out_path, threads)
    status, btindex = cleaner.buildIndex(ref_path)
    if status != 0:
        logger.error('Exiting assembler')
        sys.exit()
    #Remove contamination
    status, cread_one, cread_two, cread_sam = cleaner.deconIllumina(read_one, read_two, 
                                                                btindex)
    if status != 0:
        logger.error('Exiting assembler') 
        sys.exit()

    #Trim reads
    status, read_one, read_two = cleaner.trimIllumina(cread_one, cread_two,
                                                    adapters)
    if status != 0:
        logger.error('Exiting assembler')
        sys.exit()

    #Clean your directory young man
    cleaner.cleanSam(cread_sam)

