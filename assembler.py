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

class Assemble:
    def __init__(self, read1, read2, outdir=None, name='sample', threads=None ):
        self.read1 = os.path.abspath(read1)
        self.read2 = os.path.abspath(read2)
        if outdir != None:
            self.outdir = os.path.abspath(outdir)
        else:
            self.outdir = os.path.abspath(os.getcwd())
        if threads != None:
            self.threads = threads
        else:
            self.threads = 2
        self.name =name
        return

    def abyss(self, abyss_path='abyss-pe', abyss_param=None):
        #Create ABySS folders
        if abyss_param == None:
            abyss_param = ['k=63']
        aoutdir = '{0}/abyss'.format(self.outdir)
        if not os.path.exists(aoutdir):
            os.mkdir(aoutdir)
        #Log file path
        outlog = open('{0}/abyss.log'.format(aoutdir), 'w')
        #Prepare abyss run command
        aoutpath = 'name={0}/{1}'.format(aoutdir, self.name)
        inpath = 'in=\'{0} {1}\''.format(self.read1, self.read2)
        kparam = 'k={0}'.format(kmer)
        acmd = [abyss_path, aoutpath, inpath] + abyss_param
        print(' '.join(acmd))
        #run_prog = subprocess.Popen(cmd, stdout=outlog,
        #    stderr=outlog, shell=False)
        #Capture stdout and stderr
        #run_status = run_prog.communicate()
        run_prog.returncode = 0
        outlog.flush()
        return('{0}/{1}-contigs.fa'.format(aoutdir, self.name), run_prog.returncode)

    def ngopt(self, ngopt_path='a5_pipeline.pl', ngopt_param=None):
        if ngopt_param == None:
            ngopt_param = []
        #Create NGOPT folders
        noutdir = '{0}/ngopt'.format(self.outdir)
        if not os.path.exists(noutdir):
            os.mkdir(noutdir)
        outlog = open('{0}/ngopt.log'.format(noutdir), 'w')
        #Prepare run commands
        noutpath = '{0}/{1}'.format(noutdir, self.name)
        ncmd = [ngopt_path, self.read1, self.read2, noutpath] + ngopt_param
        #run_prog = subprocess.Popen(cmd, stdout=outlog,
        #    stderr=outlog, shell=False)
        #Capture stdout and stderr
        #run_status = run_prog.communicate()
        #outlog.flush()
        run_prog.returncode = 0
        print(' '.join(ncmd))
        return('{0}.contigs.fasta'.format(noutpath), run_prog.returncode)

    def sgaPreProcess(self, outdir=None, outlog=sys.stdout, sga_path='sga'):
        #Create sga preprocessing output directory
        if not outdir:
            outdir = self.outdir
        ppoutdir = '{0}/preprocess'.format(outdir)
        if not os.path.exists(ppoutdir):
            os.mkdir(ppoutdir)
        #Prepare preprocess commands
        ppoutpath = '{0}/{1}.fastq'.format(ppoutdir, self.name)
        ppcmd = [sga_path, 'preprocess', '-p', '1', self.read1, self.read2,
            '-o', ppoutpath]
        #run_prog = subprocess.Popen(ppcmd, stdout=outlog,
        #    stderr=outlog, shell=False)
        #Capture stdout and stderr
        #run_status = run_prog.communicate()
        #outlog.flush()
        run_prog.returncode = 0
        print(' '.join(ppcmd))
        return(ppoutpath, run_prog.returncode)

    def sgaIndex(self, ppoutpath, outdir=None, outlog=sys.stdout, sga_path='sga'):
        #Create sga index output directory
        if not outdir:
            outdir = self.outdir
        ioutdir = '{0}/index'.format(outdir)
        if not os.path.exists(ioutdir):
            os.mkdir(ioutdir)
        ioutpath =  '{0}/{1}'.format(ioutdir, self.name)
        #Prepare index command
        icmd = [sga_path, 'index', '-t', self.threads, '-a', 'ropebwt', '-p',
            ioutpath, ppoutpath]
        #run_prog = subprocess.Popen(indexcmd, stdout=outlog,
        #    stderr=outlog, shell=False)
        #Capture stdout and stderr
        #run_status = run_prog.communicate()
        #outlog.flush()
        run_prog.returncode = 0
        print(' '.join(icmd))
        return(ioutpath, run_prog.returncode)

    def sgaCorrect(self, ioutpath, ppoutpath, correct=41, outdir=None, outlog=sys.stdout, sga_path='sga'):
        #Create correction directory
        if not outdir:
            outdir = self.outdir
        coutdir = '{0}/correct'.format(outdir)
        if not os.path.exists(coutdir):
            os.mkdir(coutdir)
        coutpath = '{0}/{1}.ec.fq'.format(coutdir, self.name)
        #Prepare correction command
        ccmd = [sga_path, 'correct', '-t', self.threads, '-k', str(correct), '--learn',
            '-p', ioutpath, '-o', coutpath, ppoutpath]
        #run_prog = subprocess.Popen(ccmd, stdout=outlog,
        #    stderr=outlog, shell=False)
        #Capture stdout and stderr
        #run_status = run_prog.communicate()
        #outlog.flush()
        run_prog.returncode = 0
        print(' '.join(ccmd))
        return(coutpath, run_prog.returncode )

    def sgaFilter(self, ioutpath, coutpath, outdir=None, outlog=sys.stdout, sga_path='sga'):
        #Create filter directory
        if not outdir:
            outdir = self.outdir
        foutdir = '{0}/filter'.format(outdir)
        if not os.path.exists(foutdir):
            os.mkdir(foutdir)
        foutpath = '{0}/{1}.filter.pass.fq'.format(foutdir, self.name)
        #Prepare filter command
        fcmd = [sga_path, 'filter', '-x', '2', '-t', self.threads, '-p', ioutpath, '-o', foutpath]
        #run_prog = subprocess.Popen(fcmd, stdout=outlog,
        #    stderr=outlog, shell=False)
        #Capture stdout and stderr
        #run_status = run_prog.communicate()
        #outlog.flush()
        run_prog.returncode = 0
        print(' '.join(fcmd))
        return(foutpath, run_prog.returncode)

    def sgaOverlap(self, ioutpath, foutpath, overlap=75, outdir=None, outlog=sys.stdout, sga_path='sga'):
        #Create overlap directory
        if not outdir:
            outdir = self.outdir
        ooutdir = '{0}/overlap'.format(outdir)
        if not os.path.exists(ooutdir):
            os.mkdir(ooutdir)
        ooutpath = '{0}/{1}.ec.filter.pass.asqg.gz'.format(ooutdir, self.name)
        ocmd = [sga_path, 'overlap', '-m', str(overlap), '-t', self.threads, '-o',
                    ooutpath, '-p', ioutpath, foutpath]
        #run_prog = subprocess.Popen(ocmd, stdout=outlog,
        #    stderr=outlog, shell=False)
        #Capture stdout and stderr
        #run_status = run_prog.communicate()
        #outlog.flush()
        run_prog.returncode = 0
        print(' '.join(ocmd))
        return(ooutpath, run_prog.returncode)

    def sgaAssemble(self, ooutpath, assemble=71, outdir=None, outlog=sys.stdout, sga_path='sga'):
        #Create assemble directory
        if not outdir:
            outdir = self.outdir
        aoutdir = '{0}/assemble'.format(outdir)
        if not os.path.exists(aoutdir):
            os.mkdir(aoutdir)
        aoutpath = '{0}/{1}'.format(aoutdir, self.name)
        acmd = [sga_path, 'assemble', '-t', self.threads, '-m', str(assemble), '-o', aoutpath,
                        ooutpath]
        #run_prog = subprocess.Popen(acmd, stdout=outlog,
        #    stderr=outlog, shell=False)
        #Capture stdout and stderr
        #run_status = run_prog.communicate()
        #outlog.flush()
        run_prog.returncode = 0
        print(' '.join(acmd))
        return(aoutpath, run_prog.returncode)

    def sga(self, sga_path='sga', sga_params=None):
        #Rewrite SGA, break into parts to better allow pipelining
        if sga_params == None:
            sga_params = [41, 75, 71]
        #Create sga output directories
        soutdir = '{0}/sga'.format(self.outdir)
        if not os.path.exists(soutdir):
            os.mkdir(soutdir)
        outlog = open('{0}/sga.log'.format(soutdir), 'w')
        #Run preprocessing
        ppoutpath, ppret = Assemble.sgaPreProcess(self, soutdir, outlog, sga_path)
        #If preprocessing failed, return error code
        if ppret != 0:
            return(None, ppret)
        #Run indexing
        ioutpath, iret = Assemble.sgaIndex(self, soutdir, outlog, sga_path, ppoutpath)
        #If indexing failed, return error code
        if iret != 0:
            return(None, iret)
        #Run correction
        coutpath, cret = Assemble.sgaCorrect(self, sga_params[0], soutdir, outlog, sga_path, ioutpath, ppoutpath)
        #If correction failed, return error code
        if cret != 0:
            return(None, cret)
        #Run post correction indexing
        ioutpath, iret = Assemble.sgaIndex(self, soutdir, outlog, sga_path, coutpath)
        #If indexing fails, return error code
        if iret != 0:
            return(None, iret)
        #Run filter
        foutpath, fret = Assemble.sgaFilter(self, soutdir, outlog, sga_path, ioutpath, coutpath)
        #If filter failed, return errror code
        if fret != 0:
            return(None, fret)
        #Run overlap
        ooutpath, oret =  Assemble.sgaOverlap(self, sga_params[1], outdir, outlog, sga_path, ioutpath, foutpath)
        #If overlap failed, return error code
        if oret != 0:
            return(None, oret)
        #Run overlap
        aoutpath, aret =  Assemble.sgaOverlap(self, sga_params[2], outdir, outlog, sga_path, ioutpath, foutpath)
        #If overlap failed, return error code
        if aret != 0:
            return(None, aret)
        else:
            return('{0}-contigs.fa'.format(aoutpath), aret)

    def spades(self, spades_path='spades.py', spades_params=None):
        if spades_params == None:
            spades_params = ['--careful', '-k', '27,49,71,93,115,127']
        #Create spades directory
        soutdir = '{0}/spades'.format(self.outdir)
        if not os.path.exists(soutdir):
            os.mkdir(soutdir)
        outlog = open('{0}/spades.log'.format(soutdir), 'w')
        soutpath= '{0}/contigs.fasta'.format(soutdir)
        #Prepare spades command
        scmd = [spades_path, '--pe1-1', self.read1, '--pe1-2', self.read2,
            '-o', soutdir] + spades_params
        #run_prog = subprocess.Popen(scmd, stdout=outlog,
        #    stderr=outlog, shell=False)
        #Capture stdout and stderr
        #run_status = run_prog.communicate()
        #outlog.flush()
        run_prog.returncode = 0
        print(' '.join(scmd))
        return(soutpath, run_prog.returncode)

    def pandaseq(self, panda_path='pandaseq', panda_param=None):
        if panda_param == None:
            panda_param = ['-o', '3', '-O', '0']
        #Create pandaseq directory
        poutdir = '{0}/pandaseq'.format(self.outdir)
        if not os.path.exists(poutdir):
            os.mkdir(poutdir)
        outlog = open('{0}/pandaseq.log'.format(soutdir), 'w')
        poutpath = '{0}/{1}.fa'.format(poutdir, self.name)
        #Prepare spades command
        pcmd = [panda_path, '-f', self.read1,'-r', self.read2, '-w',
            poutpath] + panda_param
        #run_prog = subprocess.Popen(pcmd, stdout=outlog,
        #    stderr=outlog, shell=False)
        #Capture stdout and stderr
        #run_status = run_prog.communicate()
        #outlog.flush()
        run_prog.returncode = 0
        print(' '.join(pcmd))
        return(poutpath, run_prog.returncode)

def unitTest(read1, read2, name, outdir, threads):
    assemble = Assemble(read1, read2, outdir, name, threads)
    abyss_contig, returncode = assemble.abyss()
    if returncode != 0:
        print("AbySS failed")
    ngopt_contig, returncode = assemble.ngopt()
    if returncode != 0:
        print("NGOPT failed")
    sga_contig, returncode = assemble.sga()
    if returncode != 0:
        print("SGA failed")
    spades_contig, returncode = assemble.spades()
    if returncode != 0:
        print("Spades failed")
    panda_contig, returncode = assemble.pandaseq()
    if returncode != 0:
        print("Pandaseq failed")
    return

if __name__ == '__main__':
    pbrazi =  argparse.ArgumentParser(prog='assembler')
    pbrazi.add_argument('-f', '--read1', type=str, dest='read1',
        help='Input fastq files, read1 followed by read2.')
    pbrazi.add_argument('-r', '--read2', type=str, dest='read2',
        help='Input fastq files, read1 followed by read2.')
    pbrazi.add_argument('-o', '--outdir', type=str, dest='outdir', default=None,
        help='Output directory.')
    pbrazi.add_argument('-n', '--sample', type=str, dest='name',
        help='Sample name')
    pbrazi.add_argument('-t', '--threads', type=str, dest='name', default='2',
        help='Number of threads')
    pbrazi.add_argument('-m', '--mode', type=str, dest='mode',
        choices=['test'], help='Sample name')
    pbrazi.add_argument('-v', '--version', action='version', version='%(prog)s 0.9.5')
    opts = pbrazi.parse_args()
    if opts.mode == 'test':
        unitTest(opts.read1, opts.read2, opts.name, opts.outdir, opts.threads)
