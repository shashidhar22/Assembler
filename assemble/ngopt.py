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

class Ngopt:
    def __init__(self, ngopt_path, read1, read2, outdir, name, threads):
        #Initialize values and create output directories
        self.ngopt_path = ngopt_path
        self.read1 = os.path.abspath(read1)
        self.read2 = os.path.abspath(read2)
        self.outdir = '{0}/ngopt'.format(os.path.abspath(outdir))
        self.threads = threads
        self.klen = klen
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
        runlogger = open(self.runtime, 'w')
        logger = open(self.log, 'w')

        #Prepare run commands
        logger.write('Ngopt pipeline started\n')
        logger.write('Changing working directory to : {0}'.format(self.outdir))
        os.chdir(self.outdir)

        ncmd = [ngopt_path, self.read1, self.read2, noutpath] + ngopt_param
        logger.write('Running NGOPT with the following command\n')
        logger.write('{0}\n'.format(' '.join(ncmd)))

        #Run NGOPT 
        nrun = subprocess.Popen(cmd, stdout=runlogger,
            stderr=runlogger)

        #Capture stdout and stderr
        nlog = nrun.communicate()
        elapsed = timeit.default_timer() - start
        runlogger.flush()
        runlogger.flush()
        if nrun.returncode != 0 :
            logger.write('NGOPT failed with exit code : {0}; Check runtime log for details.\n'.format(nrun.returncode))
            logger.close()
            return(nrun.returncode)
        else:
            logger.write('NGOPT completed successfully; Runtime : {0}\n'.format(elapsed))
            logger.write('Contigs can be found in : {0}'.format(self.result))
            logger.close()
            return(nrun.returncode)
