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
    def __init__(self, read1, read2, outdir=os.path.abspath(os.getcwd()), name='sample'):
        self.read1 = os.path.abspath(read1)
        self.read2 = os.path.abspath(read2)
        self.outdir = os.path.abspath(outdir)
        self.name=name

    def abyss(self, kmer=63, abyss_path='abyss-pe'):
        #Create ABySS folders
        outdir = '{0}/abyss'.format(self.outdir)
        #Log file path
        outlog = open('{0}/abyss.log'.format(outdir), 'w')
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        #Prepare abyss run command
        outpath = 'name={0}/{1}'.format(outdir, self.name)
        inpath = 'in=\'{0} {1}\''.format(self.read1, self.read2)
        kparam = 'k={0}'.format(kmer)
        cmd = [abyss_path, outpath, kparam, inpath]
        run_prog = subprocess.Popen(cmd, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, shell=False)
        #Capture stdout and stderr
        run_status = run_prog.communicate()
        outlog.write('{0}\n'.format(run_status[0]))
        outlog.write('{0}\n'.format(run_status[1]))
        outlog.close()
        return('{0}-contigs.fa'.format(outpath), run_status.returncode)

    def ngopt(self, ngopt_path='a5_pipeline.pl'):
        #Create NGOPT folders
        outdir = '{0}/ngopt'.format(self.outdir)
        outlog = open('{0}/ngopt.log'.format(outdir), 'w')
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        #Prepare run commands
        outpath = '{0}/{1}'.format(outdir, self.name)
        cmd = [ngopt_path, self.read1, self.read2, outpath]
        run_prog = subprocess.Popen(cmd, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, shell=False)
        #Capture stdout and stderr
        run_status = run_prog.communicate()
        outlog.write('{0}\n'.format(run_status[0]))
        outlog.write('{1}\n'.format(run_status[1]))
        outlog.close()
        return('{0}.contigs.fasta'.format(outpath), run_status.returncode)

    def sgaPreProcess(self, threads=8, outdir=self.outdir, outlog=sys.stdout, sga_path='sga'):
        #Create sga preprocessing output directory
        ppoutdir = '{0}/preprocess'.format(outdir)
        if not os.path.exists(ppoutdir):
            os.mkdir(ppoutdir)
        #Prepare preprocess commands
        ppoutpath = '{0}/{1}.fastq'.format(ppoutdir, self.name)
        ppcmd = [sga_path, 'preprocess', '-p', '1', self.read1, self.read2,
            '-o', ppoutpath]
        run_prog = subprocess.Popen(ppcmd, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, shell=False)
        run_status = run_prog.communicate()
        outlog.write('SGA preprocessing:\n')
        outlog.write('{0}\n'.format(run_status[0]))
        outlog.write('{1}\n'.format(run_status[1]))
        return(ppoutpath, run_status.returncode)

    def sgaIndex(self, ppoutpath, threads=8, outdir=self.outdir, outlog=sys.stdout, sga_path='sga'):
        #Create sga index output directory
        ioutdir = '{0}/index'.format(outdir)
        if not os.path.exists(ioutdir):
            os.mkdir(ioutdir)
        ioutpath =  '{0}/{1}'.format(ioutdir, self.name)
        #Prepare index command
        indexcmd = [sga_path, 'index', '-t', threads, '-a', 'ropebwt', '-p',
            ioutpath, ppoutpath]
        run_prog = subprocess.Popen(indexcmd, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, shell=False)
        run_status = run_prog.communicate()
        outlog.write('{0}\n'.format(run_status[0]))
        outlog.write('{1}\n'.format(run_status[1]))
        return(ioutpath, run_status.returncode)

    def sgaCorrect(self, ioutpath, ppoutpath, threads=8, correct=41, outdir=self.outdir, outlog=sys.stdout, sga_path='sga'):
        #Create correction directory
        coutdir = '{0}/correct'.format(outdir)
        if not os.path.exists(coutdir):
            os.mkdir(coutdir)
        coutpath = '{0}/{1}.ec.fq'.format(coutdir, self.name)
        #Prepare correction command
        ccmd = [sga_path, 'correct', '-t', threads, '-k', str(correct), '--learn',
            '-p', ioutpath, '-o', coutpath, ppoutpath]
        run_prog = subprocess.Popen(ccmd, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, shell=False)
        run_status = run_prog.communicate()
        outlog.write('{0}\n'.format(run_status[0]))
        outlog.write('{1}\n'.format(run_status[1]))
        return(coutpath, run_status.returncode)

    def sgaFilter(self, ioutpath, coutpath, threads=8, outdir=self.outdir, outlog=sys.stdout, sga_path='sga'):
        #Create filter directory
        foutdir = '{0}/filter'.format(outdir)
        if not os.path.exists(foutdir):
            os.mkdir(foutdir)
        foutpath = '{0}/{1}.filter.pass.fq'.format(foutdir, self.name)
        #Prepare filter command
        fcmd = [sga_path, 'filter', '-x', '2', '-t', threads, '-p', ioutpath, '-o', foutpath]
        run_prog = subprocess.Popen(fcmd, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, shell=False)
        run_status = run_prog.communicate()
        outlog.write('{0}\n'.format(run_status[0]))
        outlog.write('{1}\n'.format(run_status[1]))
        return(foutpath, run_status.returncode)

    def sgaOverlap(self, ioutpath, foutpath, threads=8, overlap=75, outdir=self.outdir, outlog=sys.stdout, sga_path='sga'):
        #Create overlap directory
        ooutdir = '{0}/overlap'.format(outdir)
        if not os.path.exists(ooutdir):
            os.mkdir(ooutdir)
        ooutpath = '{0}/{1}.ec.filter.pass.asqg.gz'.format(ooutdir, self.name)
        ocmd = [sga_path, 'overlap', '-m', str(overlap), '-t', '8', '-o',
                    ooutpath, '-p', ioutpath, foutpath]
        run_prog = subprocess.Popen(ocmd, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, shell=False)
        run_status = run_prog.communicate()
        outlog.write('{0}\n'.format(run_status[0]))
        outlog.write('{1}\n'.format(run_status[1]))
        return(ooutpath, run_status.returncode)

    def sgaAssemble(self, ooutpath, threads=8, assemble=71, outdir=self.outdir, outlog=sys.stdout, sga_path='sga'):
        #Create assemble directory
        aoutdir = '{0}/assemble'.format(outdir)
        if not os.path.exists(aoutdir):
            os.mkdir(aoutdir)
        aoutpath = '{0}/{1}'.format(aoutdir, self.name)
        acmd = [sga_path, 'assemble', '-t', threads, '-m', str(assemble), '-o', aoutpath,
                        ooutpath]
        run_prog = subprocess.Popen(acmd, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, shell=False)
        run_status = run_prog.communicate()
        outlog.write('{0}\n'.format(run_status[0]))
        outlog.write('{1}\n'.format(run_status[1]))
        return(aoutpath, run_status.returncode)

    def sga(self, correct=41, overlap=75, assemble=71, threads=8, sga_path='sga'):
        #Rewrite SGA, break into parts to better allow pipelining
        #Create sga output directories
        outdir = '{0}/sga'.format(self.outdir)
        outlog = open('{0}/sga.log'.format(outdir), 'w')
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        #Run preprocessing
        ppoutpath, ppret = sgaPreProcess(self, threads, outdir, outlog, sga_path)
        #If preprocessing failed, return error code
        if ppret != 0:
            return(None, ppret)
        #Run indexing
        ioutpath, iret = sgaIndex(self, threads, outdir, outlog, sga_path, ppoutpath)
        #If indexing failed, return error code
        if iret != 0:
            return(None, iret)
        #Run correction
        coutpath, cret = sgaCorrect(self, threads, correct, outdir, outlog, sga_path, ioutpath, ppoutpath)
        #If correction failed, return error code
        if cret != 0:
            return(None, cret)
        #Run post correction indexing
        ioutpath, iret = sgaIndex(self, threads, outdir, outlog, sga_path, coutpath)
        #If indexing fails, return error code
        if iret != 0:
            return(None, iret)
        #Run filter
        foutpath, fret = sgaFilter(self, threads, outdir, outlog, sga_path, ioutpath, coutpath)
        #If filter failed, return errror code
        if fret != 0:
            return(None, fret)
        #Run overlap
        ooutpath, oret =  sgaOverlap(self, threads, overlap, outdir, outlog, sga_path, ioutpath, foutpath)
        #If overlap failed, return error code
        if oret != 0:
            return(None, oret)
        #Run overlap
        aoutpath, aret =  sgaOverlap(self, threads, overlap, outdir, outlog, sga_path, ioutpath, foutpath)
        #If overlap failed, return error code
        if aret != 0:
            return(None, aret)
        else:
            return('{0}-contigs.fa'.format(aoutpath), aret)

    def spades(self, k='27,49,71,93,115,127', spades_path='spades.py'):
        #Create spades directory
        soutdir = '{0}/spades'.format(self.outdir)
        if not os.path.exists(soutdir):
            os.mkdir(soutdir)
        outlog = open('{0}/spades.log'.format(soutdir), 'w')
        soutpath= '{0}/contigs.fasta'.format(soutdir)
        #Prepare spades command
        scmd = [spades_path, '--careful', '-k', ksize, '--pe1-1', self.read1,
                    '--pe1-2', self.read2, '-o', soutdir]
        run_prog = subprocess.Popen(scmd, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, shell=False)
        run_status = run_prog.communicate()
        outlog.write('{0}\n'.format(run_status[0]))
        outlog.write('{1}\n'.format(run_status[1]))
        return(soutpath, run_status.returncode)

    def pandaseq(self, minoverlap=3, maxoverlap=0, panda_path='pandaseq'):
        #Create pandaseq directory
        poutdir = '{0}/pandaseq'.format(self.outdir)
        if not os.path.exists(poutdir):
            os.mkdir(poutdir)
        outlog = open('{0}/pandaseq.log'.format(soutdir), 'w')
        poutpath = '{0}/{1}.fa'.format(poutdir, self.name)
        #Prepare spades command
        scmd = [panda_path, '-o', minoverlap, '-O', maxoverlap,
            '-f', self.read1,'-r', self.read2, '-w', poutpath]
        run_prog = subprocess.Popen(scmd, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, shell=False)
        run_status = run_prog.communicate()
        outlog.write('{0}\n'.format(run_status[0]))
        outlog.write('{1}\n'.format(run_status[1]))
        return(poutpath, run_status.returncode)

def unitTest(read1, read2, name, outdir):
    assemble = Assemble(outdir, name, read1, read2)
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
    pbrazi.add_argument('-m', '--mode', type=str, dest='mode',
        choices=['test'], help='Sample name')
    pbrazi.add_argument('-v', '--version', action='version', version='%(prog)s 0.9.5')
    opts = pbrazi.parse_args()
    if opts.mode == 'test':
        unitTest(opts.read1, opts.read2, opts.outdir, opts.name)
