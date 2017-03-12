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
from uuid import uuid4
from itertools import repeat
from multiprocessing import Pool
from collections import defaultdict
from collections import namedtuple
from assemble.readers import Fastq
from assemble.prepinputs import Prepper

logger = logging.getLogger('Nutcracker')

class Sga:
    def __init__(self, sga_path, config, sga_param, out_path, threads):
        logger = logging.getLogger('Nutcracker')
        #Initialize values and create output directories
        sga_param = sga_param.split(',')
        self.name = 'sample'
        self.sga_path = sga_path
        self.config = config
        self.out_path = '{0}/sga'.format(os.path.abspath(out_path))
        self.threads = threads
        self.preprocess = '{0}/{1}.fastq'.format(self.out_path, self.name)
        self.index = '{0}/{1}'.format(self.out_path, self.name)
        self.ckmer = sga_param[0]
        self.correct = '{0}/{1}.ec.fastq'.format(self.out_path, self.name)
        self.filter = '{0}/{1}.filter.pass.fastq'.format(self.out_path, self.name)
        self.overlap = '{0}/{1}.ec.filter.pass.asqg.gz'.format(self.out_path, self.name)
        self.okmer = sga_param[1]
        self.akmer = sga_param[2]
        self.assemble = '{0}/{1}'.format(self.out_path, self.name)
        self.results = '{0}/{1}-contigs.fa'.format(self.out_path, self.name)
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)

        return

    def prepInp(self):
        merge_fastq_rone = '{0}/{1}_r1.fastq'.format(self.out_path, self.name)
        merge_fastq_rtwo = '{0}/{1}_r2.fastq'.format(self.out_path, self.name)
        rones = open(merge_fastq_rone, 'a')
        rtwos = open(merge_fastq_rtwo, 'a')
        for samples in self.config:
            if self.config[samples].paired and self.config[samples].prep == 'Short':
                print(samples)
                frone = Fastq(self.config[samples].files[0], './', 'phred33')
                frtwo = Fastq(self.config[samples].files[1], './', 'phred33')
                for rone, rtwo in zip(frone.read(), frtwo.read()):
                    rones.write('{0}\n{1}\n{2}\n{3}\n'.format(rone.header, rone.seq, rone.sheader, rone.quals))
                    rtwos.write('{0}\n{1}\n{2}\n{3}\n'.format(rtwo.header, rtwo.seq, rtwo.sheader, rtwo.quals))
        rones.close()
        rtwos.close()
        self.rone = merge_fastq_rone
        self.rtwo = merge_fastq_rtwo
        return

    def sgaPreProcess(self):
        '''Run SGA pre processing'''
        #Start logging
        start = timeit.default_timer()

        #Prepare preprocess commands
        logger.info('SGA pipeline started')
        ppcmd = [self.sga_path, 'preprocess', '-p', '1', self.rone, self.rtwo,
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
            log = logs.decode('utf-8', errors='ignore').split('\n')
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

    def splitter(self, file_list):
        fin_list = list()
        mid = len(file_list)/2
        fin_list = [file_list[:mid], file_list[mid:]]
        return(fin_list)

    def sgaMerge(self, file_list):
        '''Run SGA Merge'''
        #Start logging
        start = timeit.default_timer()

        #Prepare correction commands
        logger.info('SGA merge started')
        logger.info('Merging {0} files'.format(len(file_list)))
        if len(file_list) == 2:
            ext = os.path.splitext(self.file_list[0])[1]
            prefix = '{0}/{1}'.format(self.out_path, uuid4(), ext)
            mcmd = [self.sga_path, 'merge', '-r', '-p', prefix, file_list[0], file_list[1]]
            logger.debug('Running SGA merge with the following command')
            logger.debug('{0}'.format(' '.join(ccmd)))

            #Run SGA correct
            mrun = subprocess.Popen(mcmd, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, shell=False)
            #Capture stdout and stderr
            mlog = mrun.communicate()
            elapsed = timeit.default_timer() - start

            for logs in mlog:
                log = logs.decode('utf-8').split('\n')
                for lines in log:
                    logger.debug(lines)

            if mrun.returncode != 0:
                logger.error('SGA merge failed with exit code : {0}; Check runtime log for details.'.format(crun.returncode))
                sys.exit()
            else:
                ret_file = glob.glob('{0}.*'.format(prefix))[0]
                logger.info('SGA merge completed successfully; Runtime : {0}'.format(elapsed))
                logger.info('Merged file can be found at : {0}'.format(ret_file))
                return([ret_file])
        elif len(file_list) > 2:
            split_list = self.splitter(file_list)
            return(sgaMerge(split_list[0]) + sgaMerge(split_list[1]))

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
        logger.info('SGA assemble started')
        acmd = [self.sga_path, 'assemble', '-t', self.threads, '-m', str(self.akmer), 
                '-o', self.assemble, self.overlap]
        logger.info('Running SGA assemble with the following command')
        logger.info('{0}'.format(' '.join(acmd)))

        #Run SGA assemble
        arun = subprocess.Popen(acmd, stdout=runlogger,
            stderr=runlogger, shell=False)
        #Capture stdout and stderr
        alog = arun.communicate()
        elapsed = timeit.default_timer() - start
        # runlogger.write('\nSGA assemble runtime logs:\n')

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
        
        #prepInp()
        preprocess_list = list()
        index_list = list()

        #Call pre process
        ret, sample_processed = self.sgaPreProcess(samples.files)
        if ret != 0:
            return(ret)

        #call index
        ret, sample_indexed = self.sgaIndex(sample_processed)
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
        ret = self.sgaOverlap()
        if ret != 0:
            return(ret)

        #call assemble
        ret = self.sgaAssemble()
        if ret != 0:
            return(ret)
        
        return(ret)

if __name__ == '__main__':
    
    #Define defaults
    sga_default = '/projects/home/sravishankar9/local/bin/sga'
    config = '/projects/home/sravishankar9/projects/Assembler/fq/test.config'
    out_path = '/projects/home/sravishankar9/projects/Assembler/local/'
    params = '41,71,75'
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
    logger.info('Running the following modes of analysis : SGA')

    sga_params = argparse.ArgumentParser(prog="SGA runner")
    sga_params.add_argument('--sga', type=str, default=sga_default,
                              help='Path to SGA executable')
    sga_params.add_argument('--config', type=str, default=config,
                              help='Path to input config file')
    sga_params.add_argument('--out_path', type=str, dest='out_path',
                              default=out_path,
                              help='Path to output directory')
    sga_params.add_argument('--params', type=str, default=params,
                              help='SGA k-mer size')
    sga_params.add_argument('--threads', type=str, default=threads,
                              help='Number of threads allocated')

    sopts = sga_params.parse_args() 
    config = main(sopts.config)

    assembler = Sga(sopts.sga, config,
                              sopts.params, sopts.out_path, sopts.threads)
    sret = assembler.sgaRun()
    contigs = assembler.results
    print(sret)
    print(contigs)

