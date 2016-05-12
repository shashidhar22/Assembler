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

class Spades:
    def __init__(self, spades_path, read1, read2, pacbio, outdir, kmers, threads ):
        #Initialize values and create output directories
        self.spades_path = spades_path
        self.read1 = os.path.abspath(read1)
        self.read2 = os.path.abspath(read2)
        self.pacbio = os.path.abspath(pacbio)
        self.outdir = '{0}/spadesHybrid'.format(os.path.abspath(outdir))
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
        runlogger = open(self.runtime, 'w')
        logger = open(self.log, 'w')

        #Prepare run commands
        logger.write('SPAdes  started\n')
        scmd = [spades_path, '--pe1-1', self.read1, '--pe1-2', self.read2,
            '--pacbio-reads', self.pacbio, '-o', self.outdir, '--careful',
            '-k', self.kmers]
        logger.write('Running SPAdes with the following command\n')
        logger.write('{0}\n'.format(' '.join(scmd)))


        #Running Spades
        srun = subprocess.Popen(srun, stdout=runlogger,
            stderr=runlogger, shell=False)
        #Capture stdout and stderr
        slog = srun.communicate()
        elapsed = timeit.default_timer() - start
        runlogger.flush()
        runlogger.close()

        if srun.returncode != 0 :
            logger.write('SPAdes hybrid failed with exit code : {0}; Check runtime log for details.\n'.format(srun.returncode))
            logger.close()
            return(srun.returncode)
        else:
            logger.write('SPAdes hybrid completed successfully; Runtime : {0}\n'.format(elapsed))
            logger.write('Contigs can be found in : {0}'.format(self.result))
            logger.close()
            return(srun.returncode)
