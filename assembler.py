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
from assemble.sprai import Sprai
from assemble.spades_pacbio import SpadesHybrid
from assemble.canu import Canu
from assemble.cleaner import Cleaner
from assemble.evaluate import Evaluate

def alignAssembly(bowtie_path, sam_path, jelly_path, read1, read2, pacbio, assembly, outdir, threads, names, xml):
    if assembly:
        for fasta, method in zip(assembly, names):
            logging.info('Aligning reads to assembly : {0}'.format(method))
            evaluate = Evaluate(bowtie_path, sam_path, jelly_path, read1, read2, None, fasta, outdir, threads, method, None)
            logging.info('Building index')
            eret = evaluate.buildIndex()
            logging.info('Aligning reads')
            eret = evaluate.alignIllumina()
            logging.info('Sorting sam file')
            eret = evaluate.sortBam()
            logging.info('Indexing bam file')
            eret = evaluate.indexBam()
            logging.info('Alignment complete for {0}'.format(method))
            evaluate.cleanSam()
    elif xml:
        for pbassembly in xml:
            logging.info('Starting PB jelly for : {0}'.format(pbassembly))
            evaluate = Evaluate(bowtie_path, sam_path, jelly_path, None, None, None, None, outdir, threads, method, pbassembly)
    return

def cleanup(bowtie_path, bbduk_path, read1, read2, pacbio, confile, outdir, threads):
    cleaner = Cleaner(bowtie_path, bbduk_path, read1, read2, pacbio, confile, outdir, threads)
    logging.info('Building bowtie index')
    cret = cleaner.buildIndex()
    if cret == 0:
        logging.info("Bowtie index built")
    else:
        logging.error("Bowtie index failed")
        sys.exit()
    
    logging.info("Removing contaminants from illumina reads")
    cret = cleaner.deconIllumina()
    if cret == 0:
        logging.info("Contaminants removed from illumina reads")
    else:
        logging.error("Bowtie failed on illumina reads")
        sys.exit()

    cleaner.cleanSam()

    logging.info("Removing contaminants from pacbio reads")
    cret = cleaner.deconPacbio()
    if cret == 0:
        logging.info("Contaminants removed from pacbio reads")
    else:
        logging.error("Bowtie failed on pacbio reads")
        sys.exit()
    cleaner.cleanSam()
    return

def realign(bowtie_path, bbduk_path, read1, read2, confile, outdir, threads):
    cleaner = Cleaner(bowtie_path, bbduk_path, read1, read2, None, confile, outdir, threads)
    logging.info('Building bowtie index')
    cret = cleaner.buildIndex()
    if cret == 0:
        logging.info("Bowtie index built")
    else:
        logging.error("Bowtie index failed")
        sys.exit()
    
    logging.info("Removing contaminants from illumina reads")
    cret = cleaner.deconIllumina()
    if cret == 0:
        logging.info("Contaminants removed from illumina reads")
    else:
        logging.error("Bowtie failed on illumina reads")
        sys.exit()
    
    logging.info("Sorting to sam files")
    cret = cleaner.sortBam()
    if cret == 0:
        logging.info("Sorting completed")
    else:
        logging.error("Sorting failed")
        sys.exit()

    logging.info("Indexing to bam files")
    cret = cleaner.indexBam()
    if cret == 0:
        logging.info("Indexing completed")
    else:
        logging.error("Indexing failed")
        sys.exit()

    cleaner.cleanSam()
    return

def unitTest(abyss_path, sga_path, spades_path, ngopt_path, panda_path, 
            blast_path, celera_path, sprai_path, canu_path, bowtie_path, 
            sam_path, bbduk_path, read1, read2, outdir, abyss_klen, 
            sga_klen, spades_klen, sample_name, threads, egs, depth, 
            reference, adapters):

    #Decontaminate illumina reads
    logging.info('Decontaminating illumina reads')
    cleaner = Cleaner(bowtie_path, sam_path, read1, read2, reference, adapters, 
                    outdir, threads) 
    logging.info('Building bowtie index')
    cleaner.buildIndex()
    logging.info('Aligning illumina reads') 
    cleaner.deconIllumina()
    logging.info('Trimming adapter and overrepresented sequences')
    cleaner.trimIllumina()
    logging.info('Decontamination of illumina reads completed')
    read1 = cleaner.tread1
    read2 = cleaner.tread2
    cleaner.cleanSam()
    logging.info('Decontaminated files can be found at : \n {0},{1}'.format(read1, read2))
    

    #AbySS call
    logging.info('Starting AbySS')
    abyss = Abyss(abyss_path, read1, read2, outdir, abyss_klen, threads)
    aret = abyss.abyss() 
    if aret != 0:
        logging.error("AbySS failed with returncode {0}; Check abyss log for further details : {1}; Exiting pipeline".format(aret, abyss.log))
        sys.exit()
    else:
        logging.info("AbySS completed sucessfully")

    #SGA call
    sga_klen = sga_klen.split(",")
    correct = sga_klen[0]
    overlap = sga_klen[1]
    assemble = sga_klen[2]
    sret = 0
    logging.info("Starting SGA")
    sga = Sga(sga_path, read1, read2, outdir, sample_name, correct, overlap, assemble, threads) 
    logging.info("Starting SGA preprocess")
    sret = sga.sgaPreProcess()
    if sret != 0:
        logging.error("SGA pre precess failed with returncode {0}; Check sga log for further details : {1}; Exiting pipeline".format(sret, sga.log))
        sys.exit()
    else:
        logging.info('SGA pre process completed sucessfully')

    logging.info("Starting SGA index")
    sret = sga.sgaIndex()
    if sret != 0:
        logging.error("SGA index failed with returncode {0}; Check sga log for further details : {1}; Exiting pipeline".format(sret, sga.log))
        sys.exit()
    else:
        logging.info("SGA index completed sucessfully")

    logging.info("Starting SGA correct")
    sret = sga.sgaCorrect()
    if sret != 0:
        logging.error("SGA correct failed with returncode {0}; Check sga log for further details : {1}; Exiting pipeline".format(sret, sga.log))
        sys.exit()
    else:
        logging.info("SGA correct completed sucessfully")
    
    logging.info("Starting SGA index")
    sret = sga.sgaIndex()
    if sret != 0:
        logging.error("SGA index failed with returncode {0}; Check sga log for further details : {1}; Exiting pipeline".format(sret, sga.log))
        sys.exit()
    else:
        logging.info("SGA index comepleted sucessfully")

    logging.info("Starting SGA filter")
    sret = sga.sgaFilter()
    if sret != 0:
        logging.error("SGA filter failed with returncode {0}; Check sga log for further details : {1}; Exiting pipeline".format(sret, sga.log))
        sys.exit()
    else:
        logging.info("SGA filter failed")

    logging.info("Starting SGA overlap")
    sret = sga.sgaOverlap()
    if sret != 0:
        logging.error("SGA overlap failed with returncode {0}; Check sga log for further details : {1}; Exiting pipeline".format(sret, sga.log))
        sys.exit()
    else:
        logging.info("SGA overlap completed sucessfully")

    logging.info("Starting SGA assemble")
    sret = sga.sgaAssemble()
    if sret != 0:
        logging.error("SGA assembly failed with returncode {0}; Check sga log for further details : {1}; Exiting pipeline".format(sret, sga.log))
        sys.exit()
    else:
        logging.info("SGA assembly comepleted sucessfully")

    logging.info("SGA pipeline completed sucessfully")

    #Spades call
    logging.info("SPAdes assembly started")
    spades = Spades(spades_path, read1, read2, outdir, spades_klen, threads)
    spret = spades.spades()
    spret = 0 
    if spret != 0:
        logging.error("SPAdes assembly failed with returncode {0}; Check spades log for further details : {1}; Exiting pipeline".format(spret, spades.log))
        sys.exit()
    else:
        logging.info("SPAdes assembly completed sucessfully")

    #Ngopt call
    logging.info("NGOPT assembly started")
    ngopt = Ngopt(ngopt_path, read1, read2, outdir, sample_name, threads)
    nret = ngopt.ngopt()
    if nret != 0:
        logging.error("NGOPT assembly failed with returncode {0}; Check NGOPT log for further details : {1}; Exiting pipeline".format(nret, ngopt.log))
        sys.exit()
    else:
        logging.info("NGOPT assembly completed sucessfully")
    
    #PandaSeq call
    logging.info("PandaSeq assembly started")   
    panda = PandaSeq(panda_path, read1, read2, outdir, threads)
    pret = panda.pandaseq()
    if pret != 0:
        logging.error("PandaSeq assembly failed with returncode {0}; Check PandaSeq log for further details : {1}; Exiting pipeline".format(pret, panda.log))
        sys.exit() 
    else:
        logging.info("PandaSeq assembly completed sucessfully")

    #Spades hybrid assembly
    logging.info("Spades hybrid assembly started")
    spades = SpadesHybrid(spades_path, read1, read2, pacbio, outdir, kmers, threads)
    spret = spades.spades()
    if spret != 0:
        logging.error("SPAdes hybrid assembly failed with returncode {0}; Check spades log for further details : {1}; Exiting pipeline".format(spret, spades.log))
        sys.exit()
    else:
        logging.info("SPAdes hybrid assembly completed sucessfully")

    #Sprai call
    logging.info("Sprai assembly started")
    sprai = Sprai(sprai_path, celera_path, blast_path, pacbio, outdir, threads, egs, depth)
    logging.info("Writing down specifications")
    sprai.ecconfig
    sprai.pacbio
    logging.info("Starting assembly")
    sret = sprai.sprai()
    if sret != 0:
        logging.error("Sprai assembly failed with returncode {0}; Check Sprai log for further details : {1}; Exiting pipeline".format(sret, sprai.log))
        sys.exit() 
    else:
        logging.info("Sprai assembly completed sucessfully")

    #Canu call
    logging.info("Canu assembly started")
    canu = Canu(canu_path, pacbio, outdir, threads, egs, depth, name)
    canu.pbconfig()
    cret = canu.canu()
    if cret != 0:
        logging.error("Canu assembly failed with returncode {0}; Check Canu log for further details : {1}; Exiting pipeline".format(cret, canu.log))
        sys.exit() 
    else:
        logging.info("Canu assembly completed sucessfully")
    
    return

if __name__ == '__main__':
    FORMAT = '%(asctime)-15s : %(levelname)-8s :  %(message)s'
    logging.basicConfig(format=FORMAT, level=logging.DEBUG)
    logging.info("Started assembly pipeline")
    pbrazi =  argparse.ArgumentParser(prog='assembler')
    pbrazi.add_argument('-f', '--read1', type=str, dest='read1', nargs='+',
        help='Input fastq files, read1 followed by read2.')
    pbrazi.add_argument('-r', '--read2', type=str, dest='read2', nargs='+',
        help='Input fastq files, read1 followed by read2.')
    pbrazi.add_argument('-p', '--pacbio', type=str, dest='pacbio', nargs='+',
        help='Input pacbio reads.')
    pbrazi.add_argument('-o', '--outdir', type=str, dest='outdir', default=None,
        help='Output directory.')
    pbrazi.add_argument('-n', '--sample', type=str, dest='name', nargs='+',
        help='Sample name')
    pbrazi.add_argument('-a', '--assembly', type=str, dest='assembly', nargs='+',
        help='Assembly file')
    pbrazi.add_argument('-t', '--threads', type=str, dest='threads', default='2',
        help='Number of threads')
    pbrazi.add_argument('-x', '--xml', type=str, dest='xml', nargs='+',
        help='Sample PB Jelly xml file')
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
    pbrazi.add_argument('--sprai', type=str, dest='sprai_path',
        help='Path to sprai', default='sprai')
    pbrazi.add_argument('--blast', type=str, dest='blast_path',
        help='Path to blastn', default='blastn')
    pbrazi.add_argument('--celera', type=str, dest='celera_path',
        help='Path to celera assembler', default='celera')
    pbrazi.add_argument('--canu', type=str, dest='canu_path',
        help='Path to canu', default='canu')
    pbrazi.add_argument('--bowtie', type=str, dest='bowtie_path',
        help='Path to bowtie', default='bowtie2')
    pbrazi.add_argument('--sam', type=str, dest='sam_path',
        help='Path to samtools', default='samtools')
    pbrazi.add_argument('--jelly', type=str, dest='jelly_path',
        help='Path to PBJelly', default='Jelly.py')
    pbrazi.add_argument('--bbduk', type=str, dest='bbduk_path',
        help='Path to bbduk', default='bbduk')
    pbrazi.add_argument('--abyss_kmers', type=str, dest='abyss_klen',
        help='Kmer length for AbySS assembly', default='63')
    pbrazi.add_argument('--sga_kmers', type=str, dest='sga_klen',
        help='Kmer length for SGA assembly', default='41,75,71')
    pbrazi.add_argument('--spaeds_kmers', type=str, dest='spades_klen', 
        help='Kmer length for Spades assembly', default='27,49,71,93,115,127')
    pbrazi.add_argument('--genome_size', type=str, dest='egs', 
        help='Estimated genome size', default='30m')
    pbrazi.add_argument('--pacbio_depth', type=int, dest='depth',
        help='Estimated depth of sequencing for pacbio read (if unknown set to 0)', default=0)
    pbrazi.add_argument('--contaminant', type=str, dest='confile',
        help='Path to contaminant reference fasta file')
    pbrazi.add_argument('--adapters', type=str, dest='adapters', nargs='+',
        help='Path to adapter and overrepresented sequence fasta file')
    pbrazi.add_argument('-m', '--mode', type=str, dest='mode',
        choices=['test', 'cleanup', 'realign', 'evaluate'], help='Sample name')
    pbrazi.add_argument('-v', '--version', action='version', version='%(prog)s 0.9.6')
    opts = pbrazi.parse_args()
    if not os.path.exists(opts.outdir):
        os.mkdir(opts.outdir)
    if opts.mod=`=jedi=3, e == 'test':=`= (abyss_path, sga_path, spades_path, ngopt_path, panda_path, blast_path, celera_path, sprai_path, canu_path, bowtie_path, sam_path, bbduk_path, read1, read2, outdir, abyss_klen, sga_klen, spades_klen, sample_name, threads, egs, depth, *_*reference*_*, adapters) =`=jedi=`=' '
        unitTest(opts.abyss_path, opts.sga_path, opts.spades_path, opts.ngopt_path, opts.panda_path,
                opts.blast_path, opts.celera_path, opts.sprai_path, opts.canu_path, 
                opts.bowtie_path, opts.sam_path, opts.bbduk_path
                opts.read1[0], opts.read2[0], opts.outdir, opts.abyss_klen, 
                opts.sga_klen, opts.spades_klen, opts.name, opts.threads, opts.egs,
                opts.depth, opts.confile, opts.adapters)

    if opts.mode == 'cleanup':
        cleanup(opts.bowtie_path, opts.bbduk_path, opts.read1[0], opts.read2[0], opts.pacbio[0], opts.confile, opts.outdir, opts.threads)
    
    if opts.mode == 'realign':
        realign(opts.bowtie_path, opts.bbduk_path, opts.read1[0], opts.read2[0], opts.confile, opts.outdir, opts.threads)

    if opts.mode == 'evaluate':
        alignAssembly(opts.bowtie_path, opts.sam_path, opts.jelly_path, opts.read1[0], opts.read2[0], None, opts.assembly, opts.outdir, opts.threads, opts.name, None)
