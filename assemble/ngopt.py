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
class Ngopt:
    def __init__(self, ngopt_path, read1, read2, outdir, name, threads):
        #Initialize values and create output directories
        self.ngopt_path = ngopt_path
        self.read1 = os.path.abspath(read1)
        self.read2 = os.path.abspath(read2)
        self.outdir = '{0}/ngopt'.format(os.path.abspath(outdir))
        self.threads = threads
        self.name = name
        self.log = '{0}/ngopt.log'.format(self.outdir)
        self.runtime = '{0}/ngopt_runtime.log'.format(self.outdir)
        self.result = '{0}/{1}.contigs.fa'.format(self.outdir, self.name)
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

        return
    
    def ngopt(self):
        '''Run NGOPT on sample'''
        #Start logging 
        start = timeit.default_timer()

        #Prepare run commands
        logger.info('Ngopt pipeline started')
        logger.debug('Changing working directory to : {0}'.format(self.outdir))
        os.chdir(self.outdir)
        runlogger = open(self.runtime, 'w')
        ncmd = [self.ngopt_path, self.read1, self.read2, './'] 
        logger.debug('Running NGOPT with the following command')
        logger.debug('{0}'.format(' '.join(ncmd)))
        #Run NGOPT 
        nrun = subprocess.Popen(ncmd, stdout=runlogger,
            stderr=runlogger)

        #Capture stdout and stderr
        nlog = nrun.communicate()
        runlogger.flush()
        runlogger.close()
        elapsed = timeit.default_timer() - start

        if nrun.returncode != 0 :
            logger.error('NGOPT failed with exit code : {0}; Check runtime log for details'.format(nrun.returncode))
            return(nrun.returncode)
        else:
            logger.info('NGOPT completed successfully; Runtime : {0}'.format(elapsed))
            logger.info('Contigs can be found in : {0}'.format(self.result))
            return(nrun.returncode)
