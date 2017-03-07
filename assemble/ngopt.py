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
    def __init__(self, ngopt_path, config, name, out_path, threads):
        #Initialize values and create output directories
        self.ngopt_path = ngopt_path

        self.rone = os.path.abspath()
        self.rtwo = os.path.abspath(rtwo)
        self.out_path = '{0}/ngopt'.format(os.path.abspath(out_path))
        self.threads = threads
        self.name = name
        self.log = '{0}/ngopt.log'.format(self.out_path)
        self.runtime = '{0}/ngopt_runtime.log'.format(self.out_path)
        self.result = '{0}/{1}.contigs.fa'.format(self.out_path, self.name)
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)

        return
    
    def ngopt(self):
        '''Run NGOPT on sample'''
        #Start logging 
        start = timeit.default_timer()

        #Prepare run commands
        logger.info('Ngopt pipeline started')
        logger.debug('Changing working directory to : {0}'.format(self.out_path))
        os.chdir(self.out_path)
        runlogger = open(self.runtime, 'w')
        ncmd = [self.ngopt_path, self.rone, self.rtwo, './'] 
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

if __name__ == '__main__':
    
    #Define defaults
    ngopt_default = '/projects/home/sravishankar9/tools/a5_miseq_linux_20150522/bin/a5_pipeline.pl'
    rone = '/projects/home/sravishankar9/projects/Assembler/fq/test_miseq_r1.fastq'
    rtwo = '/projects/home/sravishankar9/projects/Assembler/fq/test_miseq_r2.fastq'
    out_path = '/projects/home/sravishankar9/projects/Assembler/local/ngopt'
    params = 'sample'
    threads = '4'
    # create logger
    FORMAT = '%{asctime}-15s : %{levelname}-8s : %{message}s'
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
    logger.info('Running the following modes of analysis : NGOPT')

    ngopt_params = argparse.ArgumentParser(prog="NGOPT runner")
    ngopt_params.add_argument('--ngopt', type=str, default=ngopt_default,
                              help='Path to AbySS executable')
    ngopt_params.add_argument('--rone', type=str, default=rone,
                              help='Path to read one')
    ngopt_params.add_argument('--rtwo', type=str, default=rtwo,
                              help='Path to read two')
    ngopt_params.add_argument('--outdir', type=str, dest='out_path',
                              default=out_path,
                              help='Path to output directory')
    ngopt_params.add_argument('--params', type=str, default=params,
                              help='NGOPT parameters')
    ngopt_params.add_argument('--threads', type=str, default=threads,
                              help='Number of threads allocated')

    nopts = ngopt_params.parse_args() 
    assembler = Ngopt(nopts.ngopt, nopts.rone, nopts.rtwo,
                              nopts.params, nopts.out_path, nopts.threads)
    nret = assembler.ngopt()
    contigs = assembler.result
    print(nret)
    print(contigs)

