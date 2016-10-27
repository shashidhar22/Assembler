
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

        #Prepare run commands
        logger.info('Pandaseq  started\n')
        #Setup pandaseq command
        pcmd = [self.panda_path, '-f', self.read1,'-r', self.read2, '-w',
            self.result] + self.params
        logger.info('Running Pandaseq with the following command')
        logger.info('{0}'.format(' '.join(pcmd)))


        #Running pandaseq
        prun = subprocess.Popen(pcmd, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, shell=False)
    
        #Capture stdout and stderr
        plog = prun.communicate()
        elapsed = timeit.default_timer() - start
        
        for logs in plog:
            log = logs.decode('utf-8').split('\n')
            for lines in log:
                logger.debug(lines)

        if prun.returncode != 0 :
            logger.error('Pandaseq failed with exit code : {0}; Check runtime log for details'.format(prun.returncode))
            return(prun.returncode)
        else:
            logger.info('Pandaseq completed successfully; Runtime : {0}'.format(elapsed))
            logger.info('Contigs can be found in : {0}'.format(self.result))
            return(prun.returncode)
