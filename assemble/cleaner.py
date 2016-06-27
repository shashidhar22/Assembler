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
import configparser
import collections
from Bio import SeqIO
from itertools import repeat
from multiprocessing import Pool
from collections import defaultdict

class Cleaner:
    def __init__(self, bowtie_path, bbduk_path, sam_path, read1, read2, reference, adapters, outdir, threads):
        #Initialize values 
        self.bowtie_path = bowtie_path
        self.bbduk_path = bbduk_path
        self.sam_path = sam_path
        self.read1 = os.path.abspath(read1)
        self.read2 = os.path.abspath(read2)
        self.reference = os.path.abspath(reference)
        self.adapters = ','.join([os.path.abspath(vals) for vals in adapters])
        self.threads = threads
        self.outdir = '{0}/cleaned_fastq'.format(os.path.abspath(outdir))
        self.cillumina = '{0}/cleaned.fastq'.format(self.outdir)
        self.cread1 = '{0}/cleaned.1.fastq'.format(self.outdir)
        self.cread2 = '{0}/cleaned.2.fastq'.format(self.outdir)
        self.tread1 = '{0}/cleaned_trimmed_r1.fastq'.format(self.outdir)
        self.tread2 = '{0}/cleaned_trimmed_r2.fastq'.format(self.outdir)
        self.btindex  = os.path.splitext(reference)[0]
        self.log = '{0}/cleaner.log'.format(self.outdir)
        self.runtime = '{0}/cleaner_runtime.log'.format(self.outdir)
        
        #Create output directory
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)
        return
        
    def buildIndex(self):
        '''Build bowtie index if it doesnt exist'''
        #Start logging
        start = timeit.default_timer()
        runlogger = open(self.runtime, 'a')
        logger = open(self.log, 'a')
        
        #Check for index
        logger.write('Checking for index\n')
        conindex = glob.glob('{0}*bt2'.format(self.conindex))
        print(conindex)
        if conindex:
            logger.write('Bowtie index exists : {0}; Skipping step;\n'.format(self.conindex))
            logger.close()
            runlogger.close()
            return(0)
        
        #Prepare run command
        logger.write('Creating bowtie index\n')
        #Setup bowtie index command
        bicmd = ['{0}-build'.format(self.bowtie_path ), self.confile, self.conindex]
        logger.write('Running bowtie index with the following command\n')
        logger.write('{0}\n'.format(' '.join(bicmd)))
    
        #Running bowtie build
        birun = subprocess.Popen(bicmd, stdout=runlogger, stderr=runlogger)
        
        #Capture stdout and stderr
        bilog = birun.communicate()
        elapsed = timeit.default_timer() - start
        runlogger.flush()
        runlogger.close()

        if birun.returncode != 0 :
            logger.write('Bowtie index failed with exit code : {0}; Check runtime log for details.\n'.format(birun.returncode))
            logger.close()
            return(birun.returncode)
        else:
            logger.write('Bowtie index completed successfully; Runtime : {0}\n'.format(elapsed))
            logger.write('Reference index can be found at : {0}.*.bt2\n'.format(self.conindex))
            logger.close()
            return(birun.returncode)
        
        
    def deconIllumina(self):
        '''Remove contaminants from illumina reads'''
        #Start logging
        start = timeit.default_timer() 
        runlogger = open(self.runtime, 'a')
        logger = open(self.log, 'a')

        #Prepare run commands
        logger.write('Decontamination  of illumina reads  started\n')
        #Setup bowtie command
        bcmd = [self.bowtie_path, '-p', '4', '-x', self.conindex, '-1', self.read1, '-2', self.read2,
                '--un-conc', self.cillumina, '--local', '-S', '{0}.sam'.format(self.cread)]
        logger.write('Running Bowtie with the following command\n')
        logger.write('{0}\n'.format(' '.join(bcmd)))


        #Running bowtie
        brun = subprocess.Popen(bcmd, stdout=runlogger,
            stderr=runlogger, shell=False)
        
        #Capture stdout and stderr
        blog = brun.communicate()
        elapsed = timeit.default_timer() - start
        runlogger.flush()
        runlogger.close()
        
    
        if brun.returncode != 0 :
            logger.write('Bowtie failed with exit code : {0}; Check runtime log for details.\n'.format(brun.returncode))
            logger.close()
            return(brun.returncode)
        else:
            logger.write('Bowtie completed successfully; Runtime : {0}\n'.format(elapsed))
            logger.write('Decontaminated reads can be found in : {0}\n'.format(self.outdir))
            logger.close()
            return(brun.returncode)

    def trimIllumina(self):
        '''Trim adapters from illumina reads'''
        #Start logging
        start = timeit.default_timer()
        runlogger = open(self.runtime, 'a')
        logger = open(self.log, 'a')

        #Prepare run commands
        logger.write('Trimming illumina reads\n')
        #Setup bbduk command
        tcmd = [self.bbduk_path, 'ref={0}'.format(self.adapters), 'in1={0}'.format(self.cread1), 
                'in2={0}'.format(self.cread2), 'out1={0}'.format(self.tread1), 
                'out2={0}'.format(self.tread2), 'ktrim=r', 'ktrim=1', 'k=27', 'mink=11', 'qtrim=rl',
                'trimq=30', 'minlength=80']
        logger.write('Running BBDuk with the following command:\n')
        logger.write('{0}\n'.format(' '.join(tcmd)))

        #Runnig BBDuk
        trun = subprocess.Popen(tcmd, stdout=runlogger,
                stderr=runlogger, shell=False)

        #Capture the stdout and stderr
        tlog = trun.communicate()
        elapsed = timeit.default_timer() - start
        runlogger.flush()
        runlogger.close()
        
    
        if trun.returncode != 0 :
            logger.write('BBDuk failed with exit code : {0}; Check runtime log for details.\n'.format(trun.returncode))
            logger.close()
            return(trun.returncode)
        else:
            logger.write('BBDuk completed successfully; Runtime : {0}\n'.format(elapsed))
            logger.write('Decontaminated reads can be found in : {0}\n'.format(self.outdir))
            logger.close()
            return(trun.returncode)

 
    def cleanSam(self):
        '''Remove sam file'''
        os.remove('{0}.sam'.format(self.cread))
        return
