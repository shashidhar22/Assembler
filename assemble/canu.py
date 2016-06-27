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

class Canu:
    def __init__(self, canu_path,  pacbio, outdir, threads, egs, depth, name):
        #Initialize values and create output directories
        self.canu_path = canu_path
        self.name = name
        self.pacbio = os.path.abspath(pacbio)
        self.outdir = '{0}/canu'.format(os.path.abspath(outdir))
        self.threads = threads
        self.log = '{0}/canu.log'.format(self.outdir)
        self.runtime = '{0}/canu_runtime.log'.format(self.outdir)
        self.result = '{0}/{1}.contigs.fasta'.format(self.outdir, self.name)
        self.spec = '{0}/{1}.specs'.format(self.outdir, self.name)
        self.egs = egs
        self.threads = threads
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

        return

    def pbconfig(self):
        '''Write celera assembly specifications file'''
        logger = open(self.log, 'a')
        ashandle = open(self.spec, 'w')
        logger.write('Writing specs file for celera assembler\n')
        ashandle.write('''# Canu specification
        errorRate=0.01
        genomeSize={0}
        shell=/bin/sh
        java=java
        saveOverlaps=true
        saveReadCorrections=true
        saveMerCounts=true
        corOverlapper=mhap
        corReAlign=true
        corMhapSensitivity=high
        useGrid=false
        unitigger=bogart
        merylThreads={1}
        merylMemory=10
        maxThreads={1}
        maxMemory=10
        obtOverlapper=mhap
        obtReAlign=true
        obtMhapSensitivity=high
        utgOverlapper=mhap
        utgReAlign=true
        utgMhapSensitivity=high'''.format(self.egs, self.threads))
        ashandle.close()
        logger.close()

        return()

    def canu(self):
        '''Run Spades on sample'''
        #Start logging
        start = timeit.default_timer() 
        runlogger = open(self.runtime, 'w')
        logger = open(self.log, 'w')

        #Prepare run commands
        logger.write('Canu  started\n')
        ccmd = [self.canu_path, '-d', self.outdir, '-p', self.name, '-s', 
                self.spec, '-pacbio-raw', self.pacbio]
        logger.write('Running Canu with the following command\n')
        logger.write('{0}\n'.format(' '.join(scmd)))


        #Running Canu
        crun = subprocess.Popen(ccmd, stdout=runlogger,
            stderr=runlogger, shell=False)
        #Capture stdout and stderr
        clog = crun.communicate()
        elapsed = timeit.default_timer() - start
        runlogger.flush()
        runlogger.close()
        if crun.returncode != 0 :
            logger.write('Canu failed with exit code : {0}; Check runtime log for details.\n'.format(crun.returncode))
            logger.close()
            return(crun.returncode)
        else:
            logger.write('Canu completed successfully; Runtime : {0}\n'.format(elapsed))
            logger.write('Contigs can be found in : {0}'.format(self.result))
            logger.close()
            return(crun.returncode)
