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

class Assemble:
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
        

        #Setup abyss command
        abyss_param = 'k={0}'.format(self.klen)
        aoutpath = 'name={0}/sample'.format(self.outdir)
        inpath = 'in=\'{0} {1}\''.format(self.read1, self.read2)
        acmd = [abyss_path, aoutpath, inpath, abyss_param]

        #Running abyss
        arun = subprocess.Popen(' '.join(acmd), stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, shell=True)
    
        #Capture stdout and stderr
        return(run_prog.returncode)
