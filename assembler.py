import os
import re
import csv
import sys
import Bio
import glob
import time
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
from assemble.abyss import Abyss
from assemble.spades import Spades
from assemble.ngopt import Ngopt
from assemble.sga import Sga
from assemble.pandaseq import PandaSeq


def unitTest(abyss_path, sga_path, spades_path, ngopt_path, read1, read2, outdir, abyss_klen, sga_klen, spades_klen, sample_name, threads):

    #AbySS call
    abyss = Abyss(abyss_path, read1, read2, outdir, abyss_klen, threads)
    aret = abyss.abyss() 
    if aret != 0:
        print("AbySS failed with returncode {0}; Check abyss log for further details : {1}; Exiting pipeline".format(aret, abyss.log))
        sys.exit()

    #SGA call
    correct = sga_klen[0]
    overlap = sga_klen[1]
    assemble = sga_klen[2]
    sga = Sga(sga_path, read1, read2, outdir, sample_name, correct, overlap, assemble, threads) 
    sret = sga.sgaPreProcess()
    if sret != 0:
        print("SGA pre precess failed with returncode {0}; Check sga log for further details : {1}; Exiting pipeline".format(sret, sga.log))
        sys.exit()
    sret = sga.sgaIndex()
    if sret != 0:
        print("SGA index failed with returncode {0}; Check sga log for further details : {1}; Exiting pipeline".format(sret, sga.log))
        sys.exit()
    sret = sga.sgaCorrect()
    if sret != 0:
        print("SGA correct failed with returncode {0}; Check sga log for further details : {1}; Exiting pipeline".format(sret, sga.log))
        sys.exit()
    sret = sga.sgaIndex()
    if sret != 0:
        print("SGA index failed with returncode {0}; Check sga log for further details : {1}; Exiting pipeline".format(sret, sga.log))
        sys.exit()
    sret = sga.sgaFilter()
    if sret != 0:
        print("SGA filter failed with returncode {0}; Check sga log for further details : {1}; Exiting pipeline".format(sret, sga.log))
        sys.exit()
    sret = sga.Overlap()
    if sret != 0:
        print("SGA overlap failed with returncode {0}; Check sga log for further details : {1}; Exiting pipeline".format(sret, sga.log))
        sys.exit()
    sret = sga.Assemble()
    if sret != 0:
        print("SGA assembly failed with returncode {0}; Check sga log for further details : {1}; Exiting pipeline".format(sret, sga.log))
        sys.exit()
    
    #Spades call
    spades = Spades(spades_path, read1, read2, outdir, spades_klen, threads)
    spret = spades.spades()
    if spret != 0:
        print("Spades assembly failed with returncode {0}; Check spades log for further details : {1}; Exiting pipeline".format(spret, spades.log))
        sys.exit()

    #Ngopt call
    ngopt = Ngopt(ngopt_path, read1, read2, outdir, sample_name, threads)
    nret = ngopt.ngopt()
    if nret != 0:
        print("NGOPT assembly failed with returncode {0}; Check NGOPT log for further details : {1}; Exiting pipeline".format(nret, ngopt.log))
        sys.exit()

    #PandaSeq call
    panda = PandaSeq(panda_path, read1, read2, outdir, threads)
    pret = panda.pandaseq()
    if pret != 0:
        print("PandaSeq assembly failed with returncode {0}; Check PandaSeq log for further details : {1}; Exiting pipeline".format(pret, panda.log))
        sys.exit() 

    return

if __name__ == '__main__':
    pbrazi =  argparse.ArgumentParser(prog='assembler')
    pbrazi.add_argument('-f', '--read1', type=str, dest='read1', nargs='+',
        help='Input fastq files, read1 followed by read2.')
    pbrazi.add_argument('-r', '--read2', type=str, dest='read2', nargs='+',
        help='Input fastq files, read1 followed by read2.')
    pbrazi.add_argument('-p', '--pacbio', type=str, dest='pacbio', nargs='+',
        help='Input pacbio reads.')
    pbrazi.add_argument('-o', '--outdir', type=str, dest='outdir', default=None,
        help='Output directory.')
    pbrazi.add_argument('-n', '--sample', type=str, dest='name', default='sample',
        help='Sample name')
    pbrazi.add_argument('-t', '--threads', type=str, dest='threads', default='2',
        help='Number of threads')
    pbrazi.add_argument('--spades', type=str, dest='spades_path',
        help='Path to spades', default='spades.py')
    pbrazi.add_argument('--sga', type=str, dest='sga_path',
        help='Path to sga', default='sga')
    pbrazi.add_argument('--abyss', type=str, dest='abyss_path',
        help='Path to abyss.', default='abyss-pe')
    pbrazi.add_argument('--ngopt', type=str, dest='ngopt_path',
        help='Path to ngopt', default='ngopt')
    pbrazi.add_argument('--pandaseq', type=str, dest='panda_path',
        help='Path to pandaseq', default='pandaseq')
    pbrazi.add_argument('--abyss_kmers', type=str, dest='abyss_klen',
        help='Kmer length for AbySS assembly', default='63')
    pbrazi.add_argument('--sga_kmers', type=str, dest='sga_klen',
        help='Kmer length for SGA assembly', default='41,75,71')
    pbrazi.add_argument('--spaeds_kmers', type=str, dest='spades_klen', 
        help='Kmer length for Spades assembly', default='27,49,71,93,115,127')
    pbrazi.add_argument('-m', '--mode', type=str, dest='mode',
        choices=['test'], help='Sample name')
    pbrazi.add_argument('-v', '--version', action='version', version='%(prog)s 0.9.6')
    opts = pbrazi.parse_args()
    if not os.path.exists(opts.outdir):
        os.mkdir(opts.outdir)
    if opts.mode == 'test':
        unitTest(opts.read1[0], opts.read2[0], opts.pacbio[0], opts.name, opts.outdir, opts.threads)
