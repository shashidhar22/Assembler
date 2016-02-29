import os
import sys
import Bio
import argparse
import subprocess
from Bio import SeqIO
from itertools import repeat
from multiprocessing import Pool

def subSample(fastq, outdir, size=1000000, readlen=0):
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
    count  = 0
    for read1, read2 in zip(recs1, recs2):
        if size == 0 and len(read1.seq) >= readlen and len(read2.seq) >= readlen:
            SeqIO.write(read1, outfiler1, "fastq")
            SeqIO.write(read2, outfiler2, "fastq")
            count = -1
        elif count <= size and len(read1.seq) >= readlen and len(read2.seq) >= readlen:
            SeqIO.write(read1, outfiler1, "fastq")
            SeqIO.write(read2, outfiler2, "fastq")
            count += 1
        elif size == 0 and len(read1.seq) >= readlen and len(read2.seq) >= readlen:
            SeqIO.write(read1, outfiler1, "fastq")
            SeqIO.write(read2, outfiler2, "fastq")
            count = -1
        elif count > size:
            break
        else:
            continue
    return(outputr1, outputr2)

def runCorrect(arguments):
    sga_file = arguments[0]
    sga_index = arguments[1]
    ksize = arguments[2]
    outdir = arguments[3]
    outfile = '{0}/subsample_corrected.{1}.ec.fq'.format(outdir, ksize)
    metric = '{0}/subsample.{1}.metrics.txt'.format(outdir, ksize)
    sga_correct = ['sga', 'correct', '-p', sga_index, '-k', str(ksize),
                    '--learn', sga_file, '-o', outfile, '--metrics', metric]
    print('Running command : {0}'.format(' '.join(sga_correct)))
    run_sga_correct = subprocess.Popen(sga_correct, shell=False)
    run_sga_correct.wait()
    run_sga_correct = None
    return(metric)

def sgaCorrect(fastq, outdir, force, krange=[35,37,39,41,43,45], size=1000000, readlen=0):
    read1 = fastq[0]
    read2 = fastq[1]
    outdir  = os.path.abspath(outdir)
    sga_file = '{0}/sgaout.fastq'.format(outdir)
    index_file = '{0}/sgaout'.format(outdir)
    if force != True and os.path.exists(sga_file):
        print('Skipping pre preprocess')
    else:
        sga_pre_process = ['sga', 'preprocess', '-p', '1', read1, read2,
                            '>', sga_file]
        print('Running command : {0}'.format(' '.join(sga_pre_process)))
        run_pre_process = subprocess.Popen(' '.join(sga_pre_process), shell=True)
        run_pre_process.wait()
        run_pre_process = None
    if force != True and os.path.exists(index_file+'.bwt'):
        print('Skipping index')
    else:
        sga_index = ['sga', 'index', '-a', 'ropebwt', '-t', '8', '-p', index_file,
                    sga_file]
        print('Running command : {0}'.format(' '.join(sga_index)))
        run_sga_index = subprocess.Popen(sga_index,  shell=False)
        run_sga_index.wait()
        run_sga_index = None
    if force != True and os.path.exists('{0}/subsample_r1.fastq'.format(outdir)):
        print('Skipping subsampling')
    else:
        read1, read2 = subSample(fastq, outdir, size, readlen)
    sga_file = '{0}/subsample.fastq'.format(outdir)
    if force != True and os.path.exists(sga_file):
        print('Skipping subssample preprocessing')
    else:
        sga_pre_process = ['sga', 'preprocess', '-p', '1', read1, read2,
                                '>', sga_file]
        print('Running command : {0}'.format(' '.join(sga_pre_process)))
        run_pre_process = subprocess.Popen(' '.join(sga_pre_process), shell=True)
        run_pre_process.wait()
        run_pre_process = None
    pool = Pool(processes=2)
    results = pool.map(runCorrect, zip(repeat(sga_file), repeat(index_file),
                krange, repeat(outdir)))
    return

def runKan(arguments):
    fastq = [arguments[0], arguments[1]]
    outdir = arguments[2]
    ksize = arguments[3]
    outfile = '{0}/sample.{1}.kc'.format(outdir, ksize)
    tdir = '{0}/sample_{1}'.format(outdir, ksize)
    kanalyze = '/projects/home/sravishankar9/tools/kanalyze-1.0.0.dev2/count'
    kanalyze_cmd = [kanalyze, '-k', str(ksize), '--minsize', '15', '-m', 'dec',
                    '--temploc', tdir, '-o', outfile] + fastq
    print('Running Kanalyze: {0}'.format(' '.join(kanalyze_cmd)))
    run_kanalyze = subprocess.Popen(kanalyze_cmd, shell=False)
    run_kanalyze.wait()
    run_kanalyze = None
    return(outfile)

def kmerOpt(fastq, outdir, krange=[35,45,55,65,75], size=10000000, readlen=0):
    outdir = os.path.abspath(outdir)
    read1, read2 = subSample(fastq, outdir, size, readlen)
    kanalyze = '/project/home/sravishankar9/tools/kanalyze-1.0.0.dev2/count'
    pool = Pool(processes=1)
    results = pool.map(runKan, zip(repeat(read1), repeat(read2), repeat(outdir), krange))
    return

def percAnalysis(arguments):
    sam_size = arguments[0]
    read1 = arguments[1]
    read2 = arguments[2]
    outdir = arguments[3]
    ksize = arguments[4]
    outdir = '{0}/subsample{1}'.format(outdir, sam_size)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    spade_path = '/projects/home/sravishankar9/tools/SPAdes-3.7.0-Linux/bin/spades.py'
    read1, read2 = subSample([read1, read2], outdir, size=sam_size)
    spade_cmd = ['python', spade_path, '--careful', '-k', ksize, '--pe1-1', read1,
            '--pe1-2', read2, '-o', outdir]
    run_spades = subprocess.Popen(spade_cmd, shell=False)
    run_spades.wait()
    run_spades = None
    return

def multiAssemble(fastq, outdir):
    outdir = os.path.abspath(outdir)
    read1 = fastq[0]
    read2 = fastq[1]
    read_handle = SeqIO.parse(open(read1), 'fastq')
    num_reads = 0
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    for line in read_handle:
        num_reads += 1
    size_range = [int(0.1*num_reads), int(0.25*num_reads), int(0.5*num_reads),
                int(0.75*num_reads), int(0.90*num_reads), int(num_reads)]
    ksize = '27,49,71,93,115,127'
    pool = Pool(processes=2)
    subsample = pool.map(percAnalysis, zip(size_range, repeat(read1),
            repeat(read2), repeat(outdir), repeat(ksize)))
    return




if __name__ == '__main__':
    pbrazi =  argparse.ArgumentParser(prog='assembler')
    pbrazi.add_argument('-f', '--fastq', type=str, dest='fastq', nargs='+',
        help='Input fastq files, read1 followed by read2.')
    pbrazi.add_argument('-o', '--outdir', type=str, dest='outdir', default=None,
        help='Output directory.')
    pbrazi.add_argument('-k', '--krange', type=int, dest='krange', nargs='+',
        default=[35, 37, 39, 41, 43, 45], help='k-mer range to test correction')
    pbrazi.add_argument('-s', '--size', type=int, dest='samplesize',
        default=1000000, help='Number of reads subsampled. Set size == 0 for subsampling based on read length alone')
    pbrazi.add_argument('-l', '--readlen', type=int, dest='readlen',
        default=0, help='Minimum read length for each read pair.')
    pbrazi.add_argument('-m', '--mode', type=str, dest='mode',
        default='kan', choices=['sga','kan','spa'], help='Mode of optimization.')
    pbrazi.add_argument('--force', action='store_true', help='Force overwrite of results.')
    pbrazi.add_argument('-v', '--version', action='version', version='%(prog)s 0.9.0')
    opts = pbrazi.parse_args()
    if opts.mode == 'sga':
        metrics = sgaCorrect(opts.fastq, opts.outdir, opts.force, opts.krange,
                            opts.samplesize, opts.readlen)
    elif opts.modde = 'spa':
        multiAssemble(opts.fastq, opts.outdir)
    else:
        metrics = kmerOpt(opts.fastq, opts.outdir, opts.krange, opts.samplesize,
                        opts.readlen)
