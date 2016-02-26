import os
import sys
import Bio
import argparse
import subprocess
from Bio import SeqIO
from itertools import repeat
from multiprocessing import Pool

def subSample(fastq, outdir, size=1000000):
    fastqr1 = open(fastq[0])
    fastqr2 = open(fastq[1])
    recs1 = SeqIO.parse(fastqr1, 'fastq')
    recs2 = SeqIO.parse(fastqr2, 'fastq')
    outdir = os.path.abspath(outdir)
    outputr1 = '{0}/subsample_r1.fastq'.format(outdir)
    outputr2 = '{0}/subsample_r2.fastq'.format(outdir)
    outfiler1 = open(outputr1, 'w')
    outfiler2 = open(outputr2, 'w')
    print('Subsampling {0} reads from raw data'.format(size))
    for count, read1, read2 in enumerate(zip(recs1, recs2)):
        if count <= size:
            SeqIO.write(recs1, outfiler1, "fastq")
            SeqIO.write(recs2, outfiler2, "fastq")
        else:
            break
    return(outfiler1, outfiler2)

def runCorrect(sga_file, sga_index, ksize, outdir):
    outfile = '{0}/subsample_corrected.{1}.ec.fq'.format(outdir, ksize)
    metric = '{0}/subsample.{1}.metrics.txt'.format(outdir, ksize)
    sga_correct = ['sga', 'correct', '-p', sga_index, '-k', ksize,
                    '--learn', sga_file]
    print('Running command : {0}'.format(' '.join(sga_correct)))
    run_sga_correct = subprocess.Popen(sga_correct, shell=False)
    run_sga_correct.wait()
    run_sga_correct = None
    return(metric)

def sgaCorrect(fastq, outdir, krange=[35,37,39,41,43,45], size=1000000):
    read1 = fastq[0]
    read2 = fastq[1]
    outdir  = os.path.abspath(outdir)
    sga_file = '{0}/sgaout.fastq'.format(outdir)
    sga_pre_process = ['sga', 'preprocess', '-p', '1', read1, read2,
                        '>', sga_file]
    print('Running command : {0}'.format(' '.join(sga_pre_process)))
    run_pre_process = subprocess.Popen(' '.join(sga_pre_process), shell=True)
    run_pre_process.wait()
    run_pre_process = None
    sga_index = ['sga', 'index', '-a', 'ropebwt', '-t', '8',
                sga_file]
    print('Running command : {0}'.format(' '.join(sga_index)))
    run_sga_index = subprocess.Popen(sga_index,  shell=False)
    run_sga_index.wait()
    run_sga_index = None
    read1, read2 = subSample(fastq, outdir)
    sga_index = '{0}/sgaout'.format(outdir)
    sga_file = '{0}/subsample.fastq'.format(outdir)
    sga_pre_process = ['sga', 'preprocess', '-p', '1', read1, read2,
                        '>', sga_file]
    print('Running command : {0}'.format(' '.join(sga_pre_process)))
    run_pre_process = subprocess.Popen(' '.join(sga_pre_process), shell=True)
    run_pre_process.wait()
    run_pre_process = None
    pool = Pool(processes=int(len(krange)))
    results = pool.map(runCorrect, zip(repeat(sga_file), repeat(sga_index),
                krange, repeat(outdir)))

if __name__ == '__main__':
    pbrazi =  argparse.ArgumentParser(prog='assembler')
    pbrazi.add_argument('-f', '--fastq', type=str, dest='fastq', nargs='+',
        help='Input fastq files, read1 followed by read2')
    pbrazi.add_argument('-o', '--outdir', type=str, dest='outdir', default=None,
        help='Output directory')
    pbrazi.add_argument('-k', '--krange', type=int, dest='krange', nargs='+',
        default=[35, 37, 39, 41, 43, 45], help='k-mer range to test correction')
    pbrazi.add_argument('-s', '--size', type=int, dest='samplesize',
        default=1000000, help='Number of reads subsampled.')
    pbrazi.add_argument('-v', '--version', action='version', version='%(prog)s 0.9.0')
    opts = pbrazi.parse_args()
    metrics = sgaCorrect(opts.fastq, opts.outdir, opts.krange, opts.samplesize)
