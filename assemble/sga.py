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
class Sga:
    def __init__(self, sga_path, read1, read2, outdir,  sga_param, threads):
        logger = logging.getLogger('Assembler')
        #Initialize values and create output directories
        sga_param = sga_param.split(',')
        self.sga_path = sga_path
        self.read1 = os.path.abspath(read1)
        self.read2 = os.path.abspath(read2)
        self.outdir = '{0}/sga'.format(os.path.abspath(outdir))
        self.threads = threads
        self.name = sga_param[0]
        self.preprocess = '{0}/{1}.fastq'.format(self.outdir, self.name)
        self.index = '{0}/{1}'.format(self.outdir, self.name)
        self.ckmer = sga_param[1]
        self.correct = '{0}/{1}.ec.fastq'.format(self.outdir, self.name)
        self.filter = '{0}/{1}.filter.pass.fastq'.format(self.outdir, self.name)
        self.overlap = '{0}/{1}.ec.filter.pass.asqg.gz'.format(self.outdir, self.name)
        self.okmer = sga_param[2]
        self.akmer = sga_param[3]
        self.assemble = '{0}/{1}'.format(self.outdir, self.name)
        self.results = '{0}/{1}-contigs.fa'.format(self.outdir, self.name)
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

        return

    def sgaPreProcess(self):
        '''Run SGA pre processing'''
        #Start logging
        start = timeit.default_timer()

        #Prepare preprocess commands
        logger.info('SGA pipeline started')
        ppcmd = [self.sga_path, 'preprocess', '-p', '1', self.read1, self.read2,
            '-o', self.preprocess]
        logger.debug('Running SGA preprocess with the following command')
        logger.debug('{0}'.format(' '.join(ppcmd)))
        
        #Run SGA pre processing
        pprun = subprocess.Popen(ppcmd, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, shell=False)

        #Capture stdout and stderr
        pplog = pprun.communicate()
        elapsed = timeit.default_timer() - start
        for logs in pplog:
            log = logs.decode('utf-8').split('\n')
            for lines in log:
                logger.debug(lines)

        if pprun.returncode != 0:
            logger.error('SGA preprocessing failed with exit code : {0}; Check runtime log for details.'.format(pprun.returncode))
            return(pprun.returncode)
        else:
            logger.debug('SGA preprocessing completed successfully; Runtime : {0}'.format(elapsed))
            logger.debug('Preprocessed fastq file can be found at : {0}'.format(self.preprocess))
            return(pprun.returncode)


    def sgaIndex(self):
        '''Run SGA index'''
        #Start logging
        start = timeit.default_timer()

        #Prepare index commands
        logger.info('SGA index started')
        if os.path.exists(self.correct):
            icmd = [self.sga_path, 'index', '-t', self.threads, '-p',
                self.index, self.correct]
        else:
            icmd = [self.sga_path, 'index', '-t', self.threads, '-p',
                self.index, self.preprocess]
        logger.info('Running SGA index with the following command')
        logger.info('{0}'.format(' '.join(icmd)))

        #Run SGA index
        irun = subprocess.Popen(icmd, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, shell=False)

        #Capture stdout and stderr
        ilog = irun.communicate()
        elapsed = timeit.default_timer() - start
        
        for logs in ilog:
            log = logs.decode('utf-8').split('\n')
            for lines in log:
                logger.debug(lines)

        if irun.returncode != 0:
            logger.error('SGA index failed with exit code : {0}; Check runtime log for details.'.format(irun.returncode))
            return(irun.returncode)
        else:
            logger.debug('SGA index completed successfully; Runtime : {0}'.format(elapsed))
            logger.debug('Index file can be found at : {0}'.format(self.index))
            return(irun.returncode)

    def sgaCorrect(self):
        '''Run SGA Correct'''
        #Start logging
        start = timeit.default_timer()

        #Prepare correction commands
        logger.info('SGA correct started')
        ccmd = [self.sga_path, 'correct', '-t', self.threads, '-k', str(self.ckmer), '--learn',
            '-p', self.index, '-o', self.correct, self.preprocess]
        logger.debug('Running SGA correct with the following command')
        logger.debug('{0}'.format(' '.join(ccmd)))

        #Run SGA correct
        crun = subprocess.Popen(ccmd, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, shell=False)
        #Capture stdout and stderr
        clog = crun.communicate()
        elapsed = timeit.default_timer() - start

        for logs in clog:
            log = logs.decode('utf-8').split('\n')
            for lines in log:
                logger.debug(lines)

        if crun.returncode != 0:
            logger.error('SGA correct failed with exit code : {0}; Check runtime log for details.'.format(crun.returncode))
            return(crun.returncode)
        else:
            logger.info('SGA correct completed successfully; Runtime : {0}'.format(elapsed))
            logger.info('Corrected fastq file can be found at : {0}'.format(self.correct))
            return(crun.returncode)

    def sgaFilter(self):
        '''Run SGA Filter'''
        #Start logging
        start = timeit.default_timer()

        #Prepare filter commands
        logger.info('SGA filter started')
        fcmd = [self.sga_path, 'filter', '-x', '2', '-t', self.threads, '-p', self.index,
                '-o', self.filter, self.correct]
        logger.debug('Running SGA filter with the following command')
        logger.debug('{0}'.format(' '.join(fcmd)))

        #Run SGA filter
        frun = subprocess.Popen(fcmd, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, shell=False)
        #Capture stdout and stderr
        flog = frun.communicate()
        elapsed = timeit.default_timer() - start

        for logs in flog:
            log = logs.decode('utf-8').split('\n')
            for lines in log:
                logger.debug(lines)

        if frun.returncode != 0:
            logger.error('SGA filter failed with exit code : {0}; Check runtime log for details'.format(frun.returncode))
            return(frun.returncode)
        else:
            logger.info('SGA filter completed successfully; Runtime : {0}'.format(elapsed))
            logger.info('Filtered fastq file can be found at : {0}'.format(self.filter))
            return(frun.returncode)

    def sgaOverlap(self):
        '''Run SGA Overlap'''
        #Start logging
        start = timeit.default_timer()

        #Prepare overlap commands
        logger.info('SGA overlap started')
        ocmd = [self.sga_path, 'overlap', '-m', str(self.okmer), '-t', self.threads, '-o',
                    self.overlap, '-p', self.index, self.filter]
        logger.info('Running SGA overlap with the following command')
        logger.info('{0}'.format(' '.join(ocmd)))

        #Run SGA overlap
        orun = subprocess.Popen(ocmd, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, shell=False)
        #Capture stdout and stderr
        olog = orun.communicate()
        elapsed = timeit.default_timer() - start

        for logs in olog:
            log = logs.decode('utf-8').split('\n')
            for lines in log:
                logger.debug(lines)

        if orun.returncode != 0:
            logger.error('SGA overlap failed with exit code : {0}; Check runtime log for details'.format(orun.returncode))
            return(orun.returncode)
        else:
            logger.info('SGA overlap completed successfully; Runtime : {0}'.format(elapsed))
            logger.info('Overlap graph file can be found at : {0}'.format(self.overlap))
            return(orun.returncode)

    def sgaAssemble(self):
        '''Run SGA Assemble'''
        #Start logging
        start = timeit.default_timer()

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

        for logs in alog:
            log = logs.decode('utf-8').split('\n')
            for lines in log:
                logger.debug(lines)

        if arun.returncode != 0:
            logger.error('SGA assemble failed with exit code : {0}; Check runtime log for details'.format(arun.returncode))
            return(arun.returncode)
        else:
            logger.info('SGA assemble completed successfully; Runtime : {0}'.format(elapsed))
            logger.info('Assemble contig file can be found at : {0}'.format(self.results))
            return(arun.returncode)

    def sgaRun(self):
        '''Run SGA assembly pipeline'''
        #Call pre process
        ret = self.sgaPreProcess()
        if ret != 0:
            return(ret)

        #call index
        ret = self.sgaIndex()
        if ret != 0:
            return(ret)

        #call correction
        ret = self.sgaCorrect()
        if ret != 0:
            return(ret)

        #call index
        ret = self.sgaIndex()
        if ret != 0:
            return(ret)

        #call filter
        ret = self.sgaFilter()
        if ret != 0:
            return(ret)

        #call overlap
        ret = self.sgaCorrect()
        if ret != 0:
            return(ret)

        return(ret)
