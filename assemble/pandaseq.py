
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

class PandaSeq:
    def __init__(self, panda_path, read1, read2, outdir, threads ):
        #Initialize values and create output directories
        self.panda_path = panda_path
        self.read1 = os.path.abspath(read1)
        self.read2 = os.path.abspath(read2)
        self.outdir = '{0}/pandaseq'.format(os.path.abspath(outdir))
        self.threads = threads
        self.log = '{0}/pandaseq.log'.format(self.outdir)
        self.runtime = '{0}/pandaseq_runtime.log'.format(self.outdir)
        self.result = '{0}/sample.fa'.format(self.outdir)
        self.params = ['-B', '-o', '3', '-O', '0']
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

        return

    def pandaseq(self):
        '''Run PandaSeq on sample'''
        #Start logging
        start = timeit.default_timer() 
        runlogger = open(self.runtime, 'w')
        logger = open(self.log, 'w')

        #Prepare run commands
        logger.write('Pandaseq  started\n')
        #Setup pandaseq command
        pcmd = [self.panda_path, '-f', self.read1,'-r', self.read2, '-w',
            self.result] + self.params
        logger.write('Running Pandaseq with the following command\n')
        logger.write('{0}\n'.format(' '.join(pcmd)))


        #Running pandaseq
        prun = subprocess.Popen(pcmd, stdout=runlogger,
            stderr=runlogger, shell=False)
    
        #Capture stdout and stderr
        plog = prun.communicate()
        elapsed = timeit.default_timer() - start
        runlogger.flush()
        runlogger.close()

        if prun.returncode != 0 :
            logger.write('Pandaseq failed with exit code : {0}; Check runtime log for details.\n'.format(prun.returncode))
            logger.close()
            return(prun.returncode)
        else:
            logger.write('Pandaseq completed successfully; Runtime : {0}\n'.format(elapsed))
            logger.write('Contigs can be found in : {0}'.format(self.result))
            logger.close()
            return(prun.returncode)
