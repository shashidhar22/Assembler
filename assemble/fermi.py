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
from assemble.prepinputs import main
logger = logging.getLogger('Assembler')
class Fermi:
    def __init__(self, fermi_path, config, out_path, fermi_runner, threads ):
        logger = logging.getLogger('Assembler')
        #Initialize values and create output directories
        self.fermi_path = fermi_path
        self.config = config
        self.out_path = '{0}/fermi'.format(os.path.abspath(out_path))
        self.threads = threads
        self.fermi_runner = fermi_runner
        self.log = '{0}/fermi.log'.format(self.out_path)
        self.runtime = '{0}/fermi_runtime.log'.format(self.out_path)
        self.result = '{0}/fmdef.p5.fq.gz'.format(self.out_path)
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)

        return

    def prepInp(self):
        fcmd = [self.fermi_path, '-Pe', self.fermi_runner, '-p', '{0}/fmdef'.format(self.out_path),
                '-t{0}'.format(self.threads)]

        for samples in self.config:
            if self.config[samples].paired:
                fcmd += self.config[samples].files

        fcmd += ['>', '{0}/fmdef.mak'.format(self.out_path)]
        return(fcmd)


    def fermi(self):
        '''Run Fermi on sample'''
        #Start logging
        start = timeit.default_timer() 

        #Prepare run commands
        logger.info('Fermi started')
        fcmd = self.prepInp() 
        logger.debug('Running Fermi with the following command')
        logger.debug('{0}'.format(' '.join(fcmd)))


        #Running Fermi
        frun = subprocess.Popen(' '.join(fcmd), stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, shell=True)
        #Capture stdout and stderr
        flog = frun.communicate()
        elapsed = timeit.default_timer() - start
        for logs in flog:
            log = logs.decode('utf-8').split('\n')
            for lines in log:
                logger.debug(lines)
        
        if frun.returncode != 0 :
            logger.error('Fermi failed with exit code : {0}; Check runtime log for details.'.format(frun.returncode))
            return(frun.returncode)
        else:
            logger.debug('Fermi step 1 completed successfully; Runtime : {0}'.format(elapsed))
        
        #Prepare make command
        fcmd = ['make', '-f', '{0}/fmdef.mak'.format(self.out_path), '-j', self.threads]
        logger.debug('Running fermi make with the following command')
        logger.debug('{0}'.format(' '.join(fcmd)))

        #Running Fermi make
        frun = subprocess.Popen(fcmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, shell=False)
    
        #Capture stdout and stderr
        flog = frun.communicate()
        elapsed = timeit.default_timer() - start
        for logs in flog:
            log = logs.decode('utf-8').split('\n')
            for lines in log:
                logger.debug(lines)
        
        if frun.returncode != 0 :
            logger.error('Fermi failed with exit code : {0}; Check runtime log for details.'.format(frun.returncode))
            return(frun.returncode)
        else:
            logger.debug('Fermi step 2 completed successfully; Runtime : {0}'.format(elapsed))
            logger.debug('Final contig can be found at : {0}'.format(self.result))
            return(frun.returncode)

if __name__ == '__main__':
    
    #Define defaults
    fermi_default = '/projects/home/sravishankar9/tools/fermi/run-fermi.pl'
    config = '/projects/home/sravishankar9/projects/Assembler/fq/test.config'
    out_path = '/projects/home/sravishankar9/projects/Assembler/local/fermi'
    params =  '/projects/home/sravishankar9/tools/fermi/fermi'
    threads = '4'

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
    logger.info('Running the following modes of analysis : fermi')


    fermi_params = argparse.ArgumentParser(prog="Fermi runner")
    fermi_params.add_argument('--fermi', type=str, default=fermi_default,
                              help='Path to Fermi executable')
    fermi_params.add_argument('--config', type=str, default=config,
                              help='Path to input folder')
    fermi_params.add_argument('--outdir', type=str, dest='out_path',
                              default=out_path,
                              help='Path to output directory')
    fermi_params.add_argument('--params', type=str, default=params,
                              help='Fermi execution path')
    fermi_params.add_argument('--threads', type=str, default=threads,
                              help='Number of threads allocated')

    fopts = fermi_params.parse_args() 
    config = main(fopts.config)
#    input_config = open(fopts.config)
#    config = dict()
#    for lines in input_config:
#        Sample = namedtuple('Sample', ['sample', 'library', 'files', 'prep', 'paired'])
#        lines = lines.strip().split(', ')
#        if lines[0] == 'Samples':
#            continue
#        else:
#            files = glob.glob(lines[2])
#            config[lines[1]] = Sample(lines[0], lines[1], files, lines[3], int(lines[4]))
#
    assembler = Fermi(fopts.fermi, config,
                              fopts.out_path, fopts.params, fopts.threads)
    fret = assembler.fermi()
    contigs = assembler.result
    print(fret)
    print(contigs)

