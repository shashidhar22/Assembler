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

class Abyss:
    def __init__(self, abyss_path, read1, read2, outdir, klen, threads ):
        #Initialize values and create output directories
        self.abyss_path = abyss_path
        self.read1 = os.path.abspath(read1)
        self.read2 = os.path.abspath(read2)
        self.outdir = '{0}/abyss'.format(os.path.abspath(outdir))
        self.threads = threads
        self.klen = klen
        self.log = '{0}/abyss.log'.format(self.outdir)
        self.runtime = '{0}/abyss_runtime.log'.format(self.outdir)
        self.result = '{0}/sample-contigs.fa'.format(self.outdir)
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)
        return

    def abyss(self):
        '''Run ABySS on sample'''
        #Start logging
        start = timeit.default_timer() 
        runlogger = open(self.runtime, 'w')
        logger = open(self.log, 'w')

        #Prepare run commands
        logger.write('AbySS  started\n')
        #Setup abyss command
        abyss_param = 'k={0}'.format(self.klen)
        aoutpath = 'name={0}/sample'.format(self.outdir)
        inpath = 'in=\'{0} {1}\''.format(self.read1, self.read2)
        acmd = [abyss_path, aoutpath, inpath, abyss_param]
        logger.write('Running AbySS with the following command\n')
        logger.write('{0}\n'.format(' '.join(acmd)))


        #Running abyss
        arun = subprocess.Popen(' '.join(acmd), stdout=runlogger,
            stderr=runlogger, shell=True)
        
        #Capture stdout and stderr
        alog = arun.communicate()
        elapsed = timeit.default_timer() - start
        runlogger.flush()
        runlogger.close()

        if arun.returncode != 0 :
            logger.write('AbySS failed with exit code : {0}; Check runtime log for details.\n'.format(arun.returncode))
            logger.close()
            return(arun.returncode)
        else:
            logger.write('AbySS completed successfully; Runtime : {0}\n'.format(elapsed))
            logger.write('Contigs can be found in : {0}'.format(self.result))
            logger.close()
            return(arun.returncode)

