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

class Sprai:
    def __init__(self, sprai_path, celera_path, blast_path,  pacbio, outdir, threads, egs, depth):
        #Initialize values and create output directories
        self.sprai_path = sprai_path
        self.celera_path = celera_path
        self.blast_path = blast_path
        self.pacbio = os.path.abspath(pacbio)
        self.outdir = '{0}/sprai'.format(os.path.abspath(outdir))
        self.threads = threads
        self.log = '{0}/sprai.log'.format(self.outdir)
        self.runtime = '{0}/sprai_runtime.log'.format(self.outdir)
        self.result = '{0}/results/CA/9-terminator/*.scf.fasta'.format(self.outdir)
        self.ec = '{0}/ec.specs'
        self.pbasm = '{0}/pbasm.specs'
        self.egs = egs
        self.depth = depth
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

        return

    def ecconfig(self):
        '''Write error correction specification'''
        logger = open(self.log, 'a')
        echandle = open(self.ec, 'w')
        logger.write('Writing specs file for error correction\n')
        echandle.write('#Input file for sprai pipeline\n')
        echandle.write('input_for_database {0}\n'.format(self.pacbio))
        echandle.write('#Estimated genomes size\n')
        echandle.write('estimated_genome_size {0}\n'.format(egs))
        echandle.write('#Estimated coverage\n')
        echandle.write('estimated_depth 0\n')
        echandle.write('#Celera assembler path\n')
        echandle.write('ca_path {0}\n'.format(self.celera_path))
        echandle.write('#Number of processes\n')
        echandle.write('partition 12\n')
        echandle.write('#Blast path\n')
        echandle.write('blast_path {0}\n'.format(self.blast_path))
        echandle.write('#Sprai path\n')
        echandle.write('sprai_path {0}\n'.format(self.sprai_path))
        echandle.write('#Blast parameters\n')
        echandle.write('word_size 18\n')
        echandle.write('evalue 1e-50\n')
        echandle.write('num_threads {0}\n'.format(self.threads))
        echandle.write('max_target_seqs 100\n')
        echandle.write('#Trim both ends of the alignment\n')
        echandle.write('trim 42\n')
        echandle.close()
        logger.close()
        return()

    def pbasmconfig(self):
        '''Write celera assembly specifications file'''
        logger = open(self.log, 'a')
        ashandle = open(self.pbasm, 'w')
        logger.write('Writing specs file for celera assembler\n')
        ashandle.write('''# Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.\n
        # 
        # Redistribution and use in source and binary forms, with or without
        # modification, are permitted (subject to the limitations in the
        # disclaimer below) provided that the following conditions are met:
        # 
        #  * Redistributions of source code must retain the above copyright
        #    notice, this list of conditions and the following disclaimer.
        # 
        #  * Redistributions in binary form must reproduce the above
        #    copyright notice, this list of conditions and the following
        #    disclaimer in the documentation and/or other materials provided
        #    with the distribution.
        # 
        #  * Neither the name of Pacific Biosciences nor the names of its
        #    contributors may be used to endorse or promote products derived
        #    from this software without specific prior written permission.
        # 
        # NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
        # GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
        # BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
        # WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
        # OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
        # DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
        # CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
        # SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
        # LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
        # USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
        # ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
        # OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
        # OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
        # SUCH DAMAGE.
        # 
        unitigger = bogart
        #utgErrorRate = 0.015 
        #utgErrorLimit = 4.5

        cnsErrorRate = 0.25
        cgwErrorRate = 0.25
        ovlErrorRate = 0.015

        frgMinLen = 1000
        ovlMinLen = 40
          
        merSize=14

        merylMemory = 16384
        merylThreads = 8

        ovlStoreMemory = 16384

        # grid info
        useGrid = 0
        scriptOnGrid = 0
        frgCorrOnGrid = 0
        ovlCorrOnGrid = 0

        sge = -S /bin/bash -V -q all.q
        #sge = -S /bin/bash -sync y -V -q all.q
        sgeScript = -pe threads 1
        sgeConsensus = -pe threads 1
        sgeOverlap = -pe threads 4
        sgeFragmentCorrection = -pe threads 4
        sgeOverlapCorrection = -pe threads 1

        #ovlHashBits = 22
        #ovlHashBlockLength = 46871347
        #ovlRefBlockSize =  537

        ovlHashBits = 25
        ovlThreads = 4
        ovlHashBlockLength = 50000000
        ovlRefBlockSize =  100000000

        ovlConcurrency = 6
        frgCorrThreads = 4
        frgCorrBatchSize = 100000
        ovlCorrBatchSize = 100000

        cnsMinFrags = 7500
        cnsConcurrency = 24

        # change sgeName every time if you do not want to wait for the jobs not necessary to wait
        sgeName = iroha''')
        ashandle.close()
        logger.close()

        return()

    def sprai(self):
        '''Run Sprai on sample'''
        #Start logging
        start = timeit.default_timer() 
        runlogger = open(self.runtime, 'w')
        logger = open(self.log, 'w')

        #Prepare run commands
        logger.write('Srai  started\n')
        pipe_path = '{0}/ezez_vx1.pl'.format(self.sprai_path)
        scmd = [pipe_path, self.ec ,self.pbasm]
        logger.write('Running Sprai with the following command\n')
        logger.write('{0}\n'.format(' '.join(scmd)))


        #Running Spades
        srun = subprocess.Popen(scmd, stdout=runlogger,
            stderr=runlogger, shell=False)
        #Capture stdout and stderr
        slog = srun.communicate()
        elapsed = timeit.default_timer() - start
        runlogger.flush()
        runlogger.close()
        self.result = glob.glob('{0}/result*/CA/9-terminator/*.scf.fasta'.format(self.outdir))[0]
        if srun.returncode != 0 :
            logger.write('Sprai failed with exit code : {0}; Check runtime log for details.\n'.format(srun.returncode))
            logger.close()
            return(srun.returncode)
        else:
            logger.write('Sprai completed successfully; Runtime : {0}\n'.format(elapsed))
            logger.write('Contigs can be found in : {0}'.format(result))
            logger.close()
            return(srun.returncode)
