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
    def __init__(self, spades_path, read1, read2, outdir, kmers, threads ):
        logger = logging.getLogger('Assembler')
        #Initialize values and create output directories
        self.spades_path = spades_path
        self.read1 = os.path.abspath(read1)
        self.read2 = os.path.abspath(read2)
        self.outdir = '{0}/spades'.format(os.path.abspath(outdir))
        self.threads = threads
        self.kmers = kmers
        self.log = '{0}/spades.log'.format(self.outdir)
        self.runtime = '{0}/spades_runtime.log'.format(self.outdir)
        self.result = '{0}/scaffolds.fasta'.format(self.outdir)
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

        return

    def spades(self):
        '''Run Spades on sample'''
        #Start logging
        start = timeit.default_timer() 

        #Prepare run commands
        logger.info('SPAdes  started\n')
        scmd = [self.spades_path, '--pe1-1', self.read1, '--pe1-2', self.read2,
            '-o', self.outdir, '--careful', '-k', self.kmers]
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

    assembler_path = sys.argv[1]
    read_one = sys.argv[2]
    read_two = sys.argv[3]
    out_path = sys.argv[4]
    assembler_param = sys.argv[5]
    threads = '4'

    assembler = Spades(assembler_path, read_one, read_two, out_path, assembler_param, threads)
    retcode = assembler.spades()
    contigs = assembler.result
    print(contigs)

