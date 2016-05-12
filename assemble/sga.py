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

class Sga:
    def __init__(self, sga_path, read1, read2, outdir, name, correct, overlap, assemble, threads):
        #Initialize values and create output directories
        self.sga_path = sga_path
        self.read1 = os.path.abspath(read1)
        self.read2 = os.path.abspath(read2)
        self.outdir = '{0}/sga'.format(os.path.abspath(outdir))
        self.threads = threads
        self.name = name
        self.log = '{0}/sga.log'.format(self.outdir)
        self.runtime = '{0}/sga_runtime.log'.format(self.outdir)
        self.preprocess = '{0}/{1}.fastq'.format(self.outdir, self.name)
        self.index = '{0}/{1}'.format(self.outdir, self.name)
        self.ckmer = correct
        self.correct = '{0}/{1}.ec.fastq'.format(self.outdir, self.name)
        self.filter = '{0}/{1}.filter.pass.fastq'.format(self.outdir, self.name)
        self.overlap = '{0}/{1}.ec.filter.pass.asqg.gz'.format(self.outdir, self.name)
        self.okmer = overlap
        self.akmer = assemble
        self.assemble = '{0}/{1}'.format(self.outdir, self.name)
        self.results = '{0}/{1}-contigs.fa'.format(self.outdir, self.name)
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

        return

    def sgaPreProcess(self):
        '''Run SGA pre processing'''
        #Start logging
        start = timeit.default_timer()
        runlogger = open(self.runtime, 'w')
        logger = open(self.log, 'w')

        #Prepare preprocess commands
        logger.write('SGA pipeline started\n')
        ppcmd = [self.sga_path, 'preprocess', '-p', '1', self.read1, self.read2,
            '-o', self.preprocess]
        logger.write('Running SGA preprocess with the following command\n')
        logger.write('{0}\n'.format(' '.join(ppcmd)))
        
        #Run SGA pre processing
        pprun = subprocess.Popen(ppcmd, stdout=runlogger,
            stderr=runlogger, shell=False)

        #Capture stdout and stderr
        pplog = pprun.communicate()
        elapsed = timeit.default_timer() - start
        runlogger.write('\nSGA preprocessing runtime logs:\n')
        runlogger.flush()
        runlogger.close()

        if pprun.returncode != 0:
            logger.write('SGA preprocessing failed with exit code : {0}; Check runtime log for details.\n'.format(pprun.returncode))
            logger.close()
            return(pprun.returncode)
        else:
            logger.write('SGA preprocessing completed successfully; Runtime : {0}\n'.format(elapsed))
            logger.write('Preprocessed fastq file can be found at : {0}'.format(self.preprocess))
            logger.close()
            return(pprun.returncode)


    def sgaIndex(self):
        '''Run SGA index'''
        #Start logging
        start = timeit.default_timer()
        runlogger = open(self.runtime, 'a')
        logger = open(self.log, 'a')

        #Prepare index commands
        logger.write('SGA index started\n')
        if os.path.exists(self.correct):
            icmd = [self.sga_path, 'index', '-t', self.threads, '-a', 'ropebwt', '-p',
                self.index, self.correct]
        else:
            icmd = [self.sga_path, 'index', '-t', self.threads, '-a', 'ropebwt', '-p',
                self.index, self.preprocess]
        logger.write('Running SGA index with the following command\n')
        logger.write('{0}\n'.format(' '.join(icmd)))

        #Run SGA index
        irun = subprocess.Popen(icmd, stdout=runlogger,
            stderr=runlogger, shell=False)

        #Capture stdout and stderr
        ilog = irun.communicate()
        elapsed = timeit.default_timer() - start
        runlogger.write('\nSGA index runtime logs:\n')
        runlogger.flush()
        runlogger.close()

        if irun.returncode != 0:
            logger.write('SGA index failed with exit code : {0}; Check runtime log for details.\n'.format(irun.returncode))
            logger.close()
            return(irun.returncode)
        else:
            logger.write('SGA index completed successfully; Runtime : {0}\n'.format(elapsed))
            logger.write('Index file can be found at : {0}'.format(self.index))
            logger.close()
            return(irun.returncode)

    def sgaCorrect(self):
        '''Run SGA Correct'''
        #Start logging
        start = timeit.default_timer()
        runlogger = open(self.runtime, 'a')
        logger = open(self.log, 'a')

        #Prepare correction commands
        logger.write('SGA correct started\n')
        ccmd = [self.sga_path, 'correct', '-t', self.threads, '-k', str(self.ckmer), '--learn',
            '-p', self.index, '-o', self.correct, self.preprocess]
        logger.write('Running SGA correct with the following command\n')
        logger.write('{0}\n'.format(' '.join(ccmd)))

        #Run SGA correct
        crun = subprocess.Popen(ccmd, stdout=runlogger,
            stderr=runlogger, shell=False)
        #Capture stdout and stderr
        clog = crun.communicate()
        elapsed = timeit.default_timer() - start
        runlogger.write('\nSGA correct runtime logs:\n')
        runlogger.flush()
        runlogger.close()

        if crun.returncode != 0:
            logger.write('SGA correct failed with exit code : {0}; Check runtime log for details.\n'.format(crun.returncode))
            logger.close()
            return(crun.returncode)
        else:
            logger.write('SGA correct completed successfully; Runtime : {0}\n'.format(elapsed))
            logger.write('Corrected fastq file can be found at : {0}'.format(self.correct))
            logger.close()
            return(crun.returncode)

    def sgaFilter(self):
        '''Run SGA Filter'''
        #Start logging
        start = timeit.default_timer()
        runlogger = open(self.runtime, 'a')
        logger = open(self.log, 'a')

        #Prepare filter commands
        logger.write('SGA filter started\n')
        fcmd = [self.sga_path, 'filter', '-x', '2', '-t', self.threads, '-p', self.index,
                '-o', self.filter, self.correct]
        logger.write('Running SGA filter with the following command\n')
        logger.write('{0}\n'.format(' '.join(fcmd)))

        #Run SGA filter
        frun = subprocess.Popen(fcmd, stdout=runlogger,
            stderr=runlogger, shell=False)
        #Capture stdout and stderr
        flog = frun.communicate()
        elapsed = timeit.default_timer() - start
        runlogger.write('\nSGA filter runtime logs:\n')
        runlogger.flush()
        runlogger.close()

        if frun.returncode != 0:
            logger.write('SGA filter failed with exit code : {0}; Check runtime log for details.\n'.format(frun.returncode))
            logger.close()
            return(frun.returncode)
        else:
            logger.write('SGA filter  completed successfully; Runtime : {0}\n'.format(elapsed))
            logger.write('Filtered fastq file can be found at : {0}'.format(self.filter))
            logger.close()
            return(frun.returncode)

    def sgaOverlap(self):
        '''Run SGA Overlap'''
        #Start logging
        start = timeit.default_timer()
        runlogger = open(self.runtime, 'a')
        logger = open(self.log, 'a')

        #Prepare overlap commands
        logger.write('SGA overlap started\n')
        ocmd = [self.sga_path, 'overlap', '-m', str(self.okmer), '-t', self.threads, '-o',
                    self.overlap, '-p', self.index, self.filter]
        logger.write('Running SGA overlap with the following command\n')
        logger.write('{0}\n'.format(' '.join(ocmd)))

        #Run SGA overlap
        orun = subprocess.Popen(ocmd, stdout=runlogger,
            stderr=runlogger, shell=False)
        #Capture stdout and stderr
        olog = orun.communicate()
        elapsed = timeit.default_timer() - start
        runlogger.write('\nSGA overlap runtime logs:\n')
        runlogger.flush()
        runlogger.close()

        if orun.returncode != 0:
            logger.write('SGA overlap failed with exit code : {0}; Check runtime log for details.\n'.format(orun.returncode))
            logger.close()
            return(orun.returncode)
        else:
            logger.write('SGA overlap completed successfully; Runtime : {0}\n'.format(elapsed))
            logger.write('Overlap graph file can be found at : {0}'.format(self.overlap))
            logger.close()
            return(orun.returncode)

    def sgaAssemble(self):
        '''Run SGA Assemble'''
        #Start logging
        start = timeit.default_timer()
        runlogger = open(self.runtime, 'a')
        logger = open(self.log, 'a')

        #Prepare assemble commands
        logger.write('SGA assemble started\n')
        acmd = [self.sga_path, 'assemble', '-t', self.threads, '-m', str(self.akmer), 
                '-o', self.assemble, self.overlap]
        logger.write('Running SGA assemble with the following command\n')
        logger.write('{0}\n'.format(' '.join(acmd)))

        #Run SGA assemble
        arun = subprocess.Popen(acmd, stdout=runlogger,
            stderr=runlogger, shell=False)
        #Capture stdout and stderr
        alog = arun.communicate()
        elapsed = timeit.default_timer() - start
        runlogger.write('\nSGA assemble runtime logs:\n')
        runlogger.flush()
        runlogger.close()

        if arun.returncode != 0:
            logger.write('SGA assemble failed with exit code : {0}; Check runtime log for details.\n'.format(arun.returncode))
            logger.close()
            return(arun.returncode)
        else:
            logger.write('SGA assemble completed successfully; Runtime : {0}\n'.format(elapsed))
            logger.write('Assemble contig file can be found at : {0}'.format(self.results))
            logger.close()
            return(arun.returncode)
