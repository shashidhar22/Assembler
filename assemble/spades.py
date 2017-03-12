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
from collections import namedtuple
from assemble.prepinputs import Prepper

logger = logging.getLogger('Nutcracker')

class Spades:

    def __init__(self, spades_path, config, kmers, out_path, threads ):
        logger = logging.getLogger('Nutcracker')
        #Initialize values and create output directories
        self.spades_path = spades_path
        self.config = config
        self.out_path = '{0}/spades'.format(os.path.abspath(out_path))
        self.threads = threads
        self.kmers = kmers
        self.log = '{0}/spades.log'.format(self.out_path)
        self.runtime = '{0}/spades_runtime.log'.format(self.out_path)
        self.result = '{0}/scaffolds.fasta'.format(self.out_path)
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)

        return

    def prepInp(self):
        scmd = [self.spades_path]
        jump = 0
        mate = 0
        single = 0
        pacbio = 0
        for samples in self.config:
            if self.config[samples].prep == 'Short' and self.config[samples].paired:
                jump += 1
                for count, files in enumerate(self.config[samples].files, start=1):
                    scmd += ['--pe{0}-{1}'.format(jump, count), files]
            elif (self.config[samples].prep == 'Single' or self.config[samples].prep == 'Long') and (not self.config[samples].paired):
                single += 1
                scmd += ['--s{0}'.format(single), self.config[samples].files[0]]
            elif self.config[samples].prep == 'Short' and self.config[samples].paired:
                mate += 1
                for count, files in enumerate(self.config[samples].files, start=1):
                    scmd += ['--mp{0}-{1}'.format(mate, count), files]

        scmd += ['-o', self.out_path, '--careful', '-k', self.kmers]
        return(scmd)

    def spades(self):
        '''Run Spades on sample'''
        #Start logging
        start = timeit.default_timer() 

        #Prepare run commands
        logger.info('SPAdes  started\n')
        scmd = self.prepInp()
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
    out_path = '/projects/home/sravishankar9/projects/Assembler/local/spades'
    params = '23,49,71,93,115,127'
    threads = '4'

    # create logger
    logger = logging.getLogger('Nutcracker')
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
    spades_params.add_argument('--input', type=str, default=inputs,
                              help='Path to Inputs config file')
    spades_params.add_argument('--outdir', type=str, dest='out_path',
                              default=out_path,
                              help='Path to output directory')
    spades_params.add_argument('--params', type=str, default=params,
                              help='Spades k-mer size')
    spades_params.add_argument('--threads', type=str, default=threads,
                              help='Number of threads allocated')

    sopts = spades_params.parse_args() 
    prepper = Prepper(sopts.input)
    config = prepper.prepInputs(sopts.input)
    assembler = Spades(sopts.spades, config, sopts.params, sopts.out_path, sopts.threads)
    sret = assembler.spades()
    contigs = assembler.result
    print(sret)
    print(contigs)

