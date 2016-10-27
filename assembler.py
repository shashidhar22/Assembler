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
from multiprocessing import Process
from collections import defaultdict
from collections import OrderedDict
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

#For version 0.9.0 removing all pacbio paramters, this feature will be re introduced for version 0.9.1.

def illumina(values):
    assembly = values[0]
    assembler_path = values[1]
    read_one = values[2]
    read_two = values[3]
    out_path = values[4]
    assembler_param = values[5]
    threads = values[6]
    
    if assembly == 'Spades':
        assembler = Spades(assembler_path, read_one, read_two, out_path, assembler_param, threads)
        retcode = assembler.spades()
        contigs = assembler.result

    elif  assembly == 'Sga':
        assembler = Sga(assembler_path, read_one, read_two, out_path, assembler_param, threads)
        retcode = assembler.sgaRun()
        contigs = assembler.result
    elif assembly == 'Ngopt':
        assembler = Ngopt(assembler_path, read_one, read_two, out_path, assembler_param, threads)
        retcode = assembler.ngopt()
        contigs = assembler.result
    elif assembly == 'PandaSeq':
        assembler = PandaSeq(assembler_path, read_one, read_two, out_path, threads)
        retcode = assembler.pandaseq()
        contigs = assembler.result
        
    elif assembly == 'Abyss':
        assembler = Abyss(assembler_path, out_path, threads)
        retcode = assembler.abyss(read_one, read_two, assembler_param)
        contigs = assembler.result

    return(assembly, retcode, contigs)





def main(bbduk_path, bowtie_path, spades_path, sga_path, ngopt_path, 
        panda_path, abyss_path, spades_param, sga_param, ngopt_param, abyss_param, 
        read_one, read_two, ref_path, threads, mode, out_path, adapters):

    #Check if inputs exists
    read_one = os.path.abspath(read_one)
    read_two  = os.path.abspath(read_two)
    
    if not os.path.exists(read_one):
        raise FileNotFoundException('Illumina fastq file not found at {0}'.format(read_one))
    if not os.path.exists(read_two):
        raise FileNotFoundException('Illumina fastq file not found at {0}'.format(read_two))

    #Reintroduce Pacbio checks in version 0.9.1
    
    #Check for presence of output directory, if not present create one
    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    # create logger
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
    logger.info('Running the following modes of analysis : {0}'.format(','.join(mode)))

    #Check for mode of operation
    if 'prepare' in mode:
        #Get abspath of paramters
        bowtie_path = os.path.abspath(bowtie_path)  
        bbduk_path = os.path.abspath(bbduk_path)
        
        #Check for correctness of path
        if not os.path.exists(bowtie_path):
            raise FileNotFoundError('Bowtie2 not found at {0}'.format(bowtie_path))
        if not os.path.exists(bbduk_path):
            raise FileNotFoundError('BBDuk not found at {0}'.format(bbduk_path))

        #Create cleanup output path
        prep_out_path = '{0}/Prepare'.format(out_path)
        if not os.path.exists(prep_out_path):
            os.mkdir(prep_out_path)

        #Start cleaner
        cleaner = Cleaner(bowtie_path, bbduk_path, prep_out_path, threads)
        status, btindex = cleaner.buildIndex(ref_path)
        if status != 0:
            logger.error('Exiting assembler')
            sys.exit()
        #Remove contamination
        status, cread_one, cread_two, cread_sam = cleaner.deconIllumina(read_one, read_two, 
                                                                    btindex)
        if status != 0:
            logger.error('Exiting assembler') 
            sys.exit()
    
        #Trim reads
        status, read_one, read_two = cleaner.trimIllumina(cread_one, cread_two,
                                                        adapters)
        if status != 0:
            logger.error('Exiting assembler')
            sys.exit()

        #Clean your directory young man
        cleaner.cleanSam(cread_sam)

    if 'illumina' in mode:
        #Get abspath of parameters
        spades_path = os.path.abspath(spades_path)
        sga_path = os.path.abspath(sga_path)
        ngopt_path = os.path.abspath(ngopt_path)
        panda_path = os.path.abspath(panda_path)
        abyss_path = os.path.abspath(abyss_path)

        #Initialize variable:
        results = OrderedDict()

        #Check for correctness of path
        if not os.path.exists(spades_path):
            raise FileNotFoundError('SPAdes not found at {0}'.format(spades_path))
        if not os.path.exists(sga_path):
            raise FileNotFoundError('SGA not found at {0}'.format(sga_path))
        if not os.path.exists(ngopt_path):
            raise FileNotFoundError('NGOPT not found at {0}'.format(ngopt_path))
        if not os.path.exists(panda_path):
            raise FileNotFoundError('PandaSeq not found at {0}'.format(panda_path))
        if not os.path.exists(abyss_path):
            raise FileNotFoundError('ABySS not found at {0}'.format(abyss_path))

        #Create cleanup output path
        ilmn_out_path = '{0}/Illumina'.format(out_path)
        if not os.path.exists(ilmn_out_path):
            os.mkdir(ilmn_out_path)

        #Create pool of assembly processes
        #Restricting to 2 parallel process, this needs to opened up
        assembly_pool = Pool(processes=2)
        assemblers = ['Spades','Sga','Ngopt','PandaSeq','Abyss']
        assembler_path = [spades_path, sga_path, ngopt_path, panda_path, abyss_path]
        assembler_params = [spades_param, sga_param, ngopt_param, None, abyss_param]
        assembly_result = assembly_pool.map(illumina, zip(assemblers, assembler_path, repeat(read_one),
                                    repeat(read_two), repeat(ilmn_out_path), assembler_params,
                                    repeat(threads)))

        #Check if assemblies ran to completions, and add them to result dictionary
        for assembly, retcode, result in assembly_result:
            if retcode:
                raise RuntimeError('Exiting Assembler')
            else:
                results[assembly] = result

        print(results)
        return 

if __name__ == '__main__':
    #Setup logger
    FORMAT = '%(asctime)-15s : %(levelname)-8s :  %(message)s'
    logger = logging.getLogger('Nutcracker')
    logger.setLevel(logging.DEBUG)
    pbrazi =  argparse.ArgumentParser(prog='Nutcracker')
    mandatory = pbrazi.add_argument_group('Basic Configuration')
    mandatory.add_argument('-f', '--read1', type=str, dest='read1', nargs='+',
        help='Input fastq files, read1 followed by read2.')
    mandatory.add_argument('-r', '--read2', type=str, dest='read2', nargs='+',
        help='Input fastq files, read1 followed by read2.')
    mandatory.add_argument('-p', '--pacbio', type=str, dest='pacbio', nargs='+',
        help='Input pacbio reads.')
    mandatory.add_argument('-o', '--outdir', type=str, dest='out_path', default=None,
        help='Output directory.')
    mandatory.add_argument('-n', '--sample', type=str, dest='name', nargs='+',
        help='Sample name')
    mandatory.add_argument('-R', '--reference', type=str, dest='ref_path', 
        help='Reference fasta file.')
    mandatory.add_argument('-t', '--threads', type=str, dest='threads', default='2',
        help='Number of threads')
    mandatory.add_argument('-m', '--mode', type=str, dest='mode', nargs='+',
        choices=['prepare', 'illumina'], help='Sample name')
    mandatory.add_argument('-v', '--version', action='version', version='%(prog)s 0.9.6')
    jelly = pbrazi.add_argument_group('PBJelly arguments')
    jelly.add_argument('--jelly', type=str, dest='jelly_path',
        help='Path to PBJelly', default='Jelly.py')
    jelly.add_argument('-a', '--assembly', type=str, dest='assembly', nargs='+',
        help='Assembly file')
    jelly.add_argument('-x', '--xml', type=str, dest='xml', nargs='+',
        help='Sample PB Jelly xml file')
    spades = pbrazi.add_argument_group('Spades arguments')
    spades.add_argument('--spades', type=str, dest='spades_path',
        help='Path to spades', default='spades.py')
    spades.add_argument('--spades_kmers', type=str, dest='spades_param', 
        help='Kmer length for Spades assembly', default='27,49,71,93,115,127')
    sga = pbrazi.add_argument_group('SGA arguments')
    sga.add_argument('--sga', type=str, dest='sga_path',
        help='Path to sga', default='sga')
    sga.add_argument('--sga_kmers', type=str, dest='sga_param',
        help='Kmer length for SGA assembly', default='41,75,71')
    abyss = pbrazi.add_argument_group('AbySS argumetns')
    abyss.add_argument('--abyss', type=str, dest='abyss_path',
        help='Path to abyss.', default='abyss-pe')
    abyss.add_argument('--abyss_kmers', type=str, dest='abyss_param',
        help='Kmer length for AbySS assembly', default='63')
    ngopt = pbrazi.add_argument_group('NGOPT arguments')
    ngopt.add_argument('--ngopt', type=str, dest='ngopt_path',
        help='Path to ngopt', default='ngopt')
    ngopt.add_argument('--ngopt_param', type=str, dest='ngopt_param',
        help='Output prefix for NGOPT', default='sample')
    panda = pbrazi.add_argument_group('PandaSeq arguments')
    panda.add_argument('--pandaseq', type=str, dest='panda_path',
        help='Path to pandaseq', default='pandaseq')
    sprai = pbrazi.add_argument_group('Sprai arguments')
    sprai.add_argument('--sprai', type=str, dest='sprai_path',
        help='Path to sprai', default='sprai')
    blast = pbrazi.add_argument_group('Blast arguments')
    blast.add_argument('--blast', type=str, dest='blast_path',
        help='Path to blastn', default='blastn')
    celera = pbrazi.add_argument_group('Celera arguments')
    celera.add_argument('--celera', type=str, dest='celera_path',
        help='Path to celera assembler', default='celera')
    canu = pbrazi.add_argument_group('Canu arguments')
    canu.add_argument('--canu', type=str, dest='canu_path',
        help='Path to canu', default='canu')
    bowtie = pbrazi.add_argument_group('Bowtie arguments')
    bowtie.add_argument('--bowtie', type=str, dest='bowtie_path',
        help='Path to bowtie', default='bowtie2')
    samtools = pbrazi.add_argument_group('Samtools arguments')
    samtools.add_argument('--sam', type=str, dest='sam_path',
        help='Path to samtools', default='samtools')
    bbduk = pbrazi.add_argument_group('BBDuk arguments')
    bbduk.add_argument('--bbduk', type=str, dest='bbduk_path',
        help='Path to bbduk', default='bbduk')
    bbduk.add_argument('--contaminant', type=str, dest='confile',
        help='Path to contaminant reference fasta file')
    bbduk.add_argument('--adapters', type=str, dest='adapters', nargs='+',
        help='Path to adapter and overrepresented sequence fasta file')
    pacbio = pbrazi.add_argument_group('PacBio assembly arguments')
    pacbio.add_argument('--genome_size', type=str, dest='egs', 
        help='Estimated genome size', default='30m')
    pacbio.add_argument('--pacbio_depth', type=int, dest='depth',
        help='Estimated depth of sequencing for pacbio read (if unknown set to 0)', default=0)
    opts = pbrazi.parse_args()

    main(opts.bbduk_path, opts.bowtie_path, opts.spades_path, opts.sga_path, opts.ngopt_path,
        opts.panda_path, opts.abyss_path, opts.spades_param, opts.sga_param, opts.ngopt_param,
        opts.abyss_param, opts.read1[0], opts.read2[0], opts.ref_path,
        opts.threads, opts.mode, opts.out_path, opts.adapters)

