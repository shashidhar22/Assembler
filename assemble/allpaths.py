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

logger = logging.getLogger('Assembler')
class Spades:
    def __init__(self, spades_path, rone, rtwo, kmers, out_path, threads ):
        logger = logging.getLogger('Assembler')
        #Initialize values and create output directories
        self.spades_path = spades_path
        self.rone = os.path.abspath(rone)
        self.rtwo = os.path.abspath(rtwo)
        self.out_path = '{0}/spades'.format(os.path.abspath(out_path))
        self.threads = threads
        self.kmers = kmers
        self.log = '{0}/spades.log'.format(self.out_path)
        self.runtime = '{0}/spades_runtime.log'.format(self.out_path)
        self.result = '{0}/scaffolds.fasta'.format(self.out_path)
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)

        return

    def spades(self):
        '''Run Spades on sample'''
        #Start logging
        start = timeit.default_timer() 

        #Prepare run commands
        logger.info('SPAdes  started\n')
        scmd = [self.spades_path, '--pe1-1', self.rone, '--pe1-2', self.rtwo,
            '-o', self.out_path, '--careful', '-k', self.kmers]
        logger.debug('Running SPAdes with the following command\n')
        logger.debug('{0}\n'.format(' '.join(scmd)))


        #Running Spades
        srun = subprocess.Popen(scmd, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, shell=False)
        #Capture stdout and stderr
        slog = srun.communicate()
        elapsed = timeit.default_timer() - start
        for logs in slog:
            log = logs.decode('utf-8').split('\n')
            for lines in log:
                logger.debug(lines)
        
        if srun.returncode != 0 :
            logger.error('SPAdes failed with exit code : {0}; Check runtime log for details.\n'.format(srun.returncode))
            return(srun.returncode)
        else:
            logger.debug('SPAdes completed successfully; Runtime : {0}\n'.format(elapsed))
            logger.debug('Contigs can be found in : {0}'.format(self.result))
            return(srun.returncode)

if __name__ == '__main__':
    
    #Define defaults
    spades_default = '/projects/home/sravishankar9/tools/SPAdes-3.9.0-Linux/bin/spades.py'
    rone = '/projects/home/sravishankar9/projects/Assembler/fq/test_miseq_r1.fastq'
    rtwo = '/projects/home/sravishankar9/projects/Assembler/fq/test_miseq_r2.fastq'
    out_path = '/projects/home/sravishankar9/projects/Assembler/local/spades'
    params = '23,49,71,93,115,127'
    threads = '4'

    # create logger
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
    logger.info('Running the following modes of analysis : spades')


    spades_params = argparse.ArgumentParser(prog="Spades runner")
    spades_params.add_argument('--spades', type=str, default=spades_default,
                              help='Path to Spades executable')
    spades_params.add_argument('--rone', type=str, default=rone,
                              help='Path to read one')
    spades_params.add_argument('--rtwo', type=str, default=rtwo,
                              help='Path to read two')
    spades_params.add_argument('--outdir', type=str, dest='out_path',
                              default=out_path,
                              help='Path to output directory')
    spades_params.add_argument('--params', type=str, default=params,
                              help='Spades k-mer size')
    spades_params.add_argument('--threads', type=str, default=threads,
                              help='Number of threads allocated')

    sopts = spades_params.parse_args() 
    assembler = Spades(sopts.spades, sopts.rone, sopts.rtwo,
                              sopts.params, sopts.out_path, sopts.threads)
    sret = assembler.spades()
    contigs = assembler.result
    print(sret)
    print(contigs)

