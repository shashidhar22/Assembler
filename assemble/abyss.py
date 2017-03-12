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
class Abyss:
    def __init__(self, abyss_path, config, klen, outdir, threads ):
        #Initialize values and create output directories
        self.abyss_path = abyss_path
        self.config = config
        self.klen = klen
        self.out_path = '{0}/abyss'.format(os.path.abspath(outdir))
        self.threads = threads
        self.log = '{0}/abyss.log'.format(self.out_path)
        self.runtime = '{0}/abyss_runtime.log'.format(self.out_path)
        self.result = '{0}/sample-contigs.fa'.format(self.out_path)
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)
        return

    def prepInp(self):
        acmd = [self.abyss_path]
        jump = 0
        mate = 0 
        single = 0
        pacbio = 0
        keys = self.config.keys()
        libs = 'lib=\''
        mp = 'mp=\''
        longs = 'long=\''
        fastq = dict()
        for samples in keys:
            if self.config[samples].prep == 'Short' and self.config[samples].paired:
                jump += 1
                libs += 'pe{0} '.format(jump)
                fastq['pe{0}'.format(jump)] = ['pe{0}'.format(jump)] + self.config[samples].files
            elif self.config[samples].prep == 'Single' and (not self.config[samples].paired):
                single += 1
                try:
                    fastq['se'] = self.config[samples].files[0]
                except KeyError:
                    fastq['se'].append(self.config[samples].files[0])
            elif self.config[samples].prep == 'Long' and (not self.config[samples].paired):
                pacbio += 1
                longs += 'long{0} '.format(pacbio)
                fastq['long{0}'.format(pacbio)] = ['long{0}'.format(pacbio)] + self.config[samples].files
            elif self.config[samples].prep == 'Mate' and self.config[samples].paired:
                mate += 1
                mp += 'mp{0} '.format(mate)
                fastq['mp{0}'.format(mate)] = ['mp{0}'.format(mate)] + self.config[samples].files

        libs += '\''
        mp += '\''
        longs += '\''                                        
        acmd += ['k={0}'.format(self.klen), 'name={0}/sample'.format(self.out_path)]
        if libs != 'lib=\'\'':
            acmd += [libs]
        if mp != 'mp=\'\'':
            acmd += [mp]
        if longs != 'long=\'\'':
            acmd += [longs]

        for samples in fastq:
            library = ''
            for count, entities in enumerate(fastq[samples]):
                if count == 0:
                    library += '{0}=\''.format(entities)
                else:
                    library += '{0} '.format(entities)
            library += '\''
            acmd += [library]
        return(acmd)



    def abyss(self):
        '''Run ABySS on sample'''
        #Start logging
        start = timeit.default_timer() 

        #Prepare run commands
        logger.info('AbySS  started')
        #Setup abyss command

        runlogger = open(self.runtime, 'w')
        acmd = self.prepInp()
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
    
    #Define defaults
    abyss_default = '/projects/home/sravishankar9/tools/abyss/bin/abyss-pe'
    inputs = '/projects/home/sravishankar9/projects/Assembler/fq/test.config'
    out_path = '/projects/home/sravishankar9/projects/Assembler/local/abyss'
    params = 63
    threads = '4'
    # create logger
    FORMAT = '%{asctime}-15s : %{levelname}-8s : %{message}s'
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
    logger.info('Running the following modes of analysis : AbySS')

    abyss_params = argparse.ArgumentParser(prog="Abyss runner")
    abyss_params.add_argument('--abyss', type=str, default=abyss_default,
                              help='Path to AbySS executable')
    abyss_params.add_argument('--input', type=str, default=inputs,
                              help='Path to input folder')
    abyss_params.add_argument('--outdir', type=str, dest='out_path',
                              default=out_path,
                              help='Path to output directory')
    abyss_params.add_argument('--params', type=int, default=params,
                              help='Abyss k-mer size')
    abyss_params.add_argument('--threads', type=str, default=threads,
                              help='Number of threads allocated')

    aopts = abyss_params.parse_args() 
    #input_config = open(aopts.input)
    prepper = Prepper(aopts.input)
    config = prepper.prepInputs()
    assembler = Abyss(aopts.abyss, config,
                              aopts.params, aopts.out_path, aopts.threads)
    aret = assembler.abyss()
    contigs = assembler.result
    print(aret)
    print(contigs)

