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
class Abyss:
    def __init__(self, abyss_path, outdir, threads ):
        #Initialize values and create output directories
        self.abyss_path = abyss_path
        self.outdir = '{0}/abyss'.format(os.path.abspath(outdir))
        self.threads = threads
        self.log = '{0}/abyss.log'.format(self.outdir)
        self.runtime = '{0}/abyss_runtime.log'.format(self.outdir)
        self.result = '{0}/sample-contigs.fa'.format(self.outdir)
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)
        return

    def abyss(self, read1, read2, klen):
        '''Run ABySS on sample'''
        #Start logging
        start = timeit.default_timer() 

        #Prepare run commands
        logger.info('AbySS  started')
        #Setup abyss command
        abyss_param = 'k={0}'.format(klen)
        aoutpath = 'name={0}/sample'.format(self.outdir)
        inpath = 'in=\'{0} {1}\''.format(read1, read2)
        runlogger = open(self.runtime, 'w')
        acmd = [self.abyss_path, aoutpath, inpath, abyss_param]
        logger.info('Running AbySS with the following command')
        logger.info('{0}'.format(' '.join(acmd)))


        #Running abyss

        arun = subprocess.Popen(' '.join(acmd), stdout=runlogger,
            stderr=runlogger, shell=True)
        
        #Capture stdout and stderr
        alog = arun.communicate()
        elapsed = timeit.default_timer() - start
        runlogger.flush()
        runlogger.close()

        if arun.returncode != 0 :
            logger.info('AbySS failed with exit code : {0}; Check runtime log for details'.format(arun.returncode))
            return(arun.returncode)
        else:
            logger.info('AbySS completed successfully; Runtime : {0}'.format(elapsed))
            logger.info('Contigs can be found in : {0}'.format(self.result))
            return(arun.returncode)


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
    logger.info('Running the following modes of analysis : AbySS')
    assembler_path = sys.argv[1]
    read_one = sys.argv[2]
    read_two = sys.argv[3]
    out_path = sys.argv[4]
    assembler_param = int(sys.argv[5])
    threads = '4'
    assembler = Abyss(assembler_path, out_path, threads)
    retcode = assembler.abyss(read_one, read_two, assembler_param)
    contigs = assembler.result
    print(contigs)
