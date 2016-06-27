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

class Evaluate:
    def __init__(self, bowtie_path, sam_path, jelly_path, read1, read2, pacbio, assembly, outdir, threads, name, xml):
        #Initialize values and create output directories
        self.bowtie_path = bowtie_path
        self.jelly_path = jelly_path
        if read1 != None and read2 != None:
            self.read1 = os.path.abspath(read1)
            self.read2 = os.path.abspath(read2)
        if pacbio != None:
            self.pacbio = os.path.abspath(pacbio)
        self.assembly = assembly
        self.sam_path = sam_path
        self.name = name
        self.xml = xml
        self.assemblyindex = os.path.splitext(os.path.abspath(assembly))[0]
        self.outdir = '{0}/{1}'.format(os.path.abspath(outdir), name)
        self.threads = threads
        self.log = '{0}/evaluate.log'.format(self.outdir)
        self.runtime = '{0}/evaluate_runtime.log'.format(self.outdir)
        self.uread1 = '{0}/unaligned_r1.fastq'.format(self.outdir)
        self.uread2 = '{0}/unaligned_r2.fastq'.format(self.outdir)
        self.uread = '{0}/{1}'.format(self.outdir, self.name)
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
        conindex = glob.glob('{0}*bt2'.format(self.assemblyindex))
        print(conindex)
        if conindex:
            logger.write('Bowtie index exists : {0}; Skipping step;\n'.format(self.assemblyindex))
            logger.close()
            runlogger.close()
            return(0)
        
        #Prepare run command
        logger.write('Creating bowtie index\n')
        #Setup bowtie index command
        bicmd = ['{0}-build'.format(self.bowtie_path ), self.assembly, self.assemblyindex]
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
            logger.write('Decontaminated reads can be found at : {0}.*.bt2\n'.format(self.assemblyindex))
            logger.close()
            return(birun.returncode)
        
        
    def alignIllumina(self):
        '''Remove contaminants from illumina reads'''
        #Start logging
        start = timeit.default_timer() 
        runlogger = open(self.runtime, 'a')
        logger = open(self.log, 'a')

        #Prepare run commands
        logger.write('Aligning illumina reads to assembly\n')
        #Setup bowtie command
        bcmd = [self.bowtie_path, '-p', '5','-x', self.assemblyindex, '-1', self.read1, '-2', self.read2,
                '--un-conc-gz', self.uread, '--local', '-S', '{0}.sam'.format(self.uread)]
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

    def sortBam(self):
        '''Sort and index file'''
        #Start logging
        start = timeit.default_timer() 
        runlogger = open(self.runtime, 'a')
        logger = open(self.log, 'a')

        #Prepare run commands
        logger.write('Sort sam file\n')
        #Setup bowtie command
        scmd = [self.sam_path, 'sort', '--reference', self.assemblyindex, 
                '-o', '{0}.bam'.format(self.uread), '{0}.sam'.format(self.uread)]
        logger.write('Running Samtools sort with the following command\n')
        logger.write('{0}\n'.format(' '.join(scmd)))


        #Running samtools
        srun = subprocess.Popen(scmd, stdout=runlogger,
            stderr=runlogger, shell=False)
        
        #Capture stdout and stderr
        slog = srun.communicate()
        elapsed = timeit.default_timer() - start
        runlogger.flush()
        runlogger.close()
        
    
        if srun.returncode != 0 :
            logger.write('Samtools sort failed with exit code : {0}; Check runtime log for details.\n'.format(srun.returncode))
            logger.close()
            return(srun.returncode)
        else:
            logger.write('Samtools sort completed successfully; Runtime : {0}\n'.format(elapsed))
            logger.write('Sorted bam can be found in : {0}\n'.format(self.outdir))
            logger.close()
            return(srun.returncode)

    def indexBam(self):
        '''Index file'''
        #Start logging
        start = timeit.default_timer() 
        runlogger = open(self.runtime, 'a')
        logger = open(self.log, 'a')

        #Prepare run commands
        logger.write('Index sam file\n')
        #Setup bowtie command
        scmd = [self.sam_path, 'index', '{0}.bam'.format(self.uread)]
        logger.write('Running Samtools index with the following command\n')
        logger.write('{0}\n'.format(' '.join(scmd)))


        #Running samtools
        srun = subprocess.Popen(scmd, stdout=runlogger,
            stderr=runlogger, shell=False)
        
        #Capture stdout and stderr
        slog = srun.communicate()
        elapsed = timeit.default_timer() - start
        runlogger.flush()
        runlogger.close()
        
    
        if srun.returncode != 0 :
            logger.write('Samtools index failed with exit code : {0}; Check runtime log for details.\n'.format(srun.returncode))
            logger.close()
            return(srun.returncode)
        else:
            logger.write('Samtools index completed successfully; Runtime : {0}\n'.format(elapsed))
            logger.write('Indexed bam can be found in : {0}\n'.format(self.outdir))
            logger.close()
            return(srun.returncode)
 
    def jellyPacbio(self):
        '''Run pb jelly on pacbio reads'''
        #Start logging
        start = timeit.default_timer() 
        runlogger = open(self.runtime, 'a')
        logger = open(self.log, 'a')

        #Prepare run commands
        logger.write('Setup pbjelly\n')
        #Setup pbjelly command
        secmd = [self.jelly_path, 'setup' ,self.xml]
        macmd = [self.jelly_path, 'mapping', self.xml]
        sucmd = [self.jelly_path, 'support', self.xml]
        excmd = [self.jelly_path, 'extraction', self.xml]
        ascmd = [self.jelly_path, 'assembly', self.xml]
        otcmd = [self.jelly_path, 'output', self.xml]
        logger.write('Running  with the following command\n')
        logger.write('{0}\n'.format(' '.join(secmd)))
        logger.write('{0}\n'.format(' '.join(macmd)))
        logger.write('{0}\n'.format(' '.join(sucmd)))
        logger.write('{0}\n'.format(' '.join(excmd)))
        logger.write('{0}\n'.format(' '.join(ascmd)))
        logger.write('{0}\n'.format(' '.join(otcmd)))
        


        #Running setup
        pbrun = subprocess.Popen(secmd, stdout=runlogger,
            stderr=runlogger, shell=False)
        pbrun.wait()
        #Running mapping
        pbrun = subprocess.Popen(macmd, stdout=runlogger,
            stderr=runlogger, shell=False)
        pbrun.wait()
        #Running support
        pbrun = subprocess.Popen(sucmd, stdout=runlogger,
            stderr=runlogger, shell=False)
        pbrun.wait()
        #Running extraction
        pbrun = subprocess.Popen(excmd, stdout=runlogger,
            stderr=runlogger, shell=False)
        pbrun.wait()
        #Running assembly
        pbrun = subprocess.Popen(ascmd, stdout=runlogger,
            stderr=runlogger, shell=False)
        pbrun.wait()
        #Running output
        pbrun = subprocess.Popen(otcmd, stdout=runlogger,
            stderr=runlogger, shell=False)
        pbrun.wait()
        return

    def cleanSam(self):
        '''Remove sam file'''
        os.remove('{0}.sam'.format(self.uread))
        return

