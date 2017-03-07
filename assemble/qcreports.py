import os
import sys
import csv
import time
import timeit
import logging
import argparse
import subprocess


class Fastqc:

    def __init__(self, fastqc, rone, rtwo, out_path, sample, threads):

        self.fastqc = os.path.abspath(fastqc)
        self.rone = os.path.abspath(rone)
        self.rtwo = os.path.abspath(rtwo)
        self.out_path = os.path.abspath(out_path)
        self.threads = str(threads)
        if sample == None or sample == ' ':
            self.sam_name = rone.split('_')[0]
        else:            
            self.sam_name =  sample
        
        return

    def runFastqc(self):
        '''Run Fastqc on the given sample'''
        
        qc_out_path = '{0}/{1}'.format(self.out_path, self.sam_name)
        if not os.path.exists(qc_out_path):
            os.mkdir(qc_out_path)
        
        qcmd = [self.fastqc, '-o', qc_out_path, '-f', 'fastq', '-t', self.threads,
                '-q', self.rone, self.rtwo]

        qrun = subprocess.Popen(qcmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, shell=False)
        
        qlogs = qrun.communicate()
        qret = qrun.returncode
        if qret != 0:
            print('Fastqc failed')
        print(' '.join(qcmd))

        return

    def run

        return


if __name__ == '__main__':

    #Setup logger
    FORMAT = '%{asctime}-15s : %{levelname}-8s : %{message}s'
    logger = logging.getLogger('QcReporter')
    logger.setLevel(logging.DEBUG)
    nqc  = argparse.ArgumentParser(prog='NutcrackerQC')
    mandatory = nqc.add_argument_group('Basic configuration')
    mandatory.add_argument('-q', '--fastqc', type=str, dest='fastqc_path',
                           help='Path to Fastqc executable.')
    mandatory.add_argument('-f', '--rone', type=str, dest='rone', help='Path \
                           to read one of paired ended fastq files.')
    mandatory.add_argument('-r', '--rtwo', type=str, dest='rtwo', help='Path \
                           to read two of paired ended fastq files.')
    mandatory.add_argument('-o', '--outdir', type=str, dest='out_path',
                           help='Path to output directory.')
    mandatory.add_argument('-n', '--sample', type=str, dest='sam_name',
                          help='Sample name. By default, the basename of \
                          read one will be used as the sample name.')
    mandatory.add_argument('-t', '--threads', type=int, dest='threads',
                           help='Number of used to run Fastqc.')
    nqc_args = nqc.parse_args()

    qcreporter = Fastqc(nqc_args.fastqc_path, nqc_args.rone, nqc_args.rtwo,
                        nqc_args.out_path, nqc_args.sam_name, nqc_args.threads)
    qcreporter.runFastqc()



    

        


        

                

        

        


