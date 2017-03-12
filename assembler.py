import os
import timeit
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
from assemble.fermi import Fermi
from assemble.prepinputs import Prepper
#For version 0.9.0 removing all pacbio paramters, this feature will be re introduced for version 0.9.1.

def illumina(values):
    assembly = values[0]
    assembler_path = values[1]
    config = values[2]
    out_path = values[3]
    assembler_param = values[4]
    threads = values[5]
    retcode = 0
    contigs = None

    if assembly == 'Spades':
        print(assembler_path, config, assembler_param, out_path, threads)
        assembler = Spades(assembler_path, config, assembler_param, out_path, threads)
        retcode = assembler.spades()
        contigs = assembler.result

    elif  assembly == 'Sga':
        print(assembler_path, config, assembler_param, out_path, threads)
        assembler = Sga(assembler_path, config, assembler_param, out_path, threads)
        retcode = assembler.sgaRun()
        contigs = assembler.results

   elif assembly == 'Ngopt':
       assembler = Ngopt(assembler_path, config, assembler_param, out_path, threads)
       retcode = assembler.ngopt()
       contigs = assembler.result

    elif assembly == 'Fermi':
        print(assembler_path, config, assembler_param, out_path, threads)
        assembler = Fermi(assembler_path, config, assembler_param, out_path, threads)
        retcode = assembler.fermi()
        contigs = assembler.result

    elif assembly == 'Abyss':
        print(assembler_path, config, assembler_param, out_path, threads)
        assembler = Abyss(assembler_path, config, assembler_param, out_path, threads)
        retcode = assembler.abyss()
        contigs = assembler.result

    return(assembly, retcode, contigs)



def prepIllumina(assembly):

    if assembly == "Spades":
        for value in range(21, 31, 4):
            params = [value + 16 for i in range(0,5)]
        yield(params)

    if assembly == 'Sga':
        #41,75,71
        for kone, ktwo, ktri in zip(range(41,51,2), range(65,75,2), range(67,77,2)):
            params = [kone, ktwo, ktri]
            yield(params)

    if assembly == 'Abyss':
        for kmer in (61,71,2)
            yield(kmer)




def main(bbduk_path, bowtie_path, spades_path, sga_path, ngopt_path, 
        fermi_path, abyss_path, spades_param, sga_param, ngopt_param, 
        fermi_param, abyss_param, input_path, ref_path, threads, mode,
        out_path, adapters, name):

    #Check if inputs exists
    config = input_path

    #Check for presence of output directory, if not present create one
    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    # create logger
    logger = logging.getLogger('Nutcracker')
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

        #Build bowtie index for host genome if it doesnt exist
        #bowtie2-build <reference path> <bowtie index path>
        status, btindex = cleaner.buildIndex(ref_path)
        if status != 0:
            logger.error('Exiting assembler')
            sys.exit()

        #Remove contamination
        #Iterate through each sample in the config dictionary and decontaminate the reads
        for samples in config:
            if config[samples].paired:
                read_one = config[samples].files[0]
                read_two = config[samples].files[1]

                #Align reads to host genome reference, and remove any read that maps to the reference genome
                #Bowtie2 run with default parameters 
                #bowtie2 -p <threads> -x <btindex> -1 <read one> -2 <read two> --un-conc <path to cleaned reads> -S <output samfile>
                status, cread_one, cread_two, cread_sam = cleaner.deconIllumina(read_one, read_two, 
                    btindex)
                if status != 0:
                    logger.error('Exiting assembler') 
                    sys.exit()
        
                #Trim reads using bbduk and list of adadpter sequences supplied by the user
                #bbduk.sh ref=<adapters.fa> in=<path to input read1> in2=<path to input read2> out=<path to cleaned read1> out2=<path to cleaned read2> ktrim=r k=27 mink=11 qtrim=rl trimq=30 minlength=80 overwrite=t
                status, read_one, read_two = cleaner.trimIllumina(cread_one, cread_two,
                                                            adapters)
                if status != 0:
                    logger.error('Exiting assembler')
                    sys.exit()

                #Delete sam files and intermedieate alignment files
                cleaner.cleanSam(cread_sam)
                config[samples].files[0] = read_one
                config[samples].files[1] = read_two

            else:
                read_one = config[samples].files
                read_two = None

                #Align reads to host genome reference, and remove any read that maps to the reference genome
                #Bowtie2 run with default parameters 
                #bowtie2 -p <threads> -x <btindex> -U < comma separated list ofsingle ended reads> --un-conc <path to cleaned reads> -S <output samfile>
                status, cread_one, cread_two, cread_sam = cleaner.deconIlluminaSinle(read_one, read_two, btindex)
                if status != 0:
                    logger.error('Exiting assembler') 
                    sys.exit()
        
                #Trim reads
                status, read_one, read_two = cleaner.trimIllumina(cread_one, cread_two,
                                                            adapters)
                if status != 0:
                    logger.error('Exiting assembler')
                    sys.exit()

                #Delete sam files and intermedieate alignment files
                cleaner.cleanSam(cread_sam)
                config[samples].files = read_one
                

    if 'illumina' in mode:
        #Get abspath of parameters
        spades_path = os.path.abspath(spades_path)
        sga_path = os.path.abspath(sga_path)
        fermi_path = os.path.abspath(fermi_path)
        abyss_path = os.path.abspath(abyss_path)

        #Initialize variable:
        results = OrderedDict()

        #Check for correctness of path
        if not os.path.exists(spades_path):
            raise FileNotFoundError('SPAdes not found at {0}'.format(spades_path))
        if not os.path.exists(sga_path):
            raise FileNotFoundError('SGA not found at {0}'.format(sga_path))
        if not os.path.exists(ngopt_path):
            raise FileNotFoundError('A5 MiSeq assembly pipeline not found at {0}'.format(ngopt_path))
        if not os.path.exists(fermi_path):
            raise FileNotFoundError('Fermi not found at {0}'.format(fermi_path))
        if not os.path.exists(abyss_path):
            raise FileNotFoundError('ABySS not found at {0}'.format(abyss_path))

        #Create cleanup output path
        ilmn_out_path = '{0}/Illumina'.format(out_path)
        if not os.path.exists(ilmn_out_path):
            os.mkdir(ilmn_out_path)

        #Create pool of assembly processes
        #Restricting to 2 parallel process, this needs to opened up
        assembly_pool = Pool(processes=2)
        assemblers = ['Spades','Sga','Fermi','Abyss']
        assembler_path = [spades_path, sga_path, fermi_path, abyss_path]

        #If optimize mode is selected, prepare kmer parameters for all assemblers tested
        if 'optimize' in mode:
            assembler_commands = list()

            for assembler, paths in zip(assemblers, assembler_path):
                assembler_optimizer = self.prepIllumina(assembler)
                for params in assembler_optimizer:
                    try:
                        assembler_commands.append([assembler, paths, config, ilmn_out_path, params, threads])
                    except ValueError:
                        assembler_commands = [[assembler, paths, config, ilmn_out_path, params, threads]]                                        

            assembly_result = assembly_pool.map(illumina, assembler_commands)

        #Else use default assembler parameters
        else:
            assembler_params = [spades_param, sga_param, fermi_param,                         abyss_param]
            assembly_result = assembly_pool.map(illumina, zip(assemblers, assembler_path, repeat(config), repeat(ilmn_out_path), assembler_params, repeat(threads)))            

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
    mandatory.add_argument('-i', '--input', type=str, dest='input_path', default=None,
        help='Path to input folder')                           
    mandatory.add_argument('-o', '--outdir', type=str, dest='out_path', default=None,
        help='Path to output directory')
    mandatory.add_argument('-n', '--sample', type=str, dest='name', nargs='+',
        help='Sample name for the project')
    mandatory.add_argument('-R', '--reference', type=str, dest='ref_path', 
        help='Path to reference fasta file')
    mandatory.add_argument('-t', '--threads', type=str, dest='threads', default='2',
        help='Number of threads for analysis')
    mandatory.add_argument('-m', '--mode', type=str, dest='mode', nargs='+',
        choices=['prepare', 'illumina', 'optimize'], help='Mode of operation. Multiple options maybe specified to create a workflow')
    mandatory.add_argument('-v', '--version', action='version', version='%(prog)s 0.9.6')
    spades = pbrazi.add_argument_group('Arguments for SPAdes analysis')
    spades.add_argument('--spades', type=str, dest='spades_path',
        help='Path to SPAdes executable', default='spades.py')
    spades.add_argument('--spades_kmers', type=str, dest='spades_param', 
        help='Kmer length for Spades assembly', default='27,49,71,93,115,127')
    sga = pbrazi.add_argument_group('Arguments for SGA analysis')
    sga.add_argument('--sga', type=str, dest='sga_path',
        help='Path to SGA executable', default='sga')
    sga.add_argument('--sga_kmers', type=str, dest='sga_param',
        help='Kmer length for SGA assembly', default='41,75,71')
    abyss = pbrazi.add_argument_group('Arguments for AbySS')
    abyss.add_argument('--abyss', type=str, dest='abyss_path',
        help='Path to ABySS executable', default='abyss-pe')
    abyss.add_argument('--abyss_kmers', type=str, dest='abyss_param',
        help='Kmer length for AbySS assembly', default='63')
    ngopt = pbrazi.add_argument_group('Arguments for A5 MiSeq assembly pipeline')
    ngopt.add_argument('--ngopt', type=str, dest='ngopt_path',
        help='Path to A5 MiSeq assembly executable', default='ngopt')
    ngopt.add_argument('--ngopt_param', type=str, dest='ngopt_param',
        help='Output prefix for A5 assembly pipeline', default='sample')
    fermi = pbrazi.add_argument_group('Arguments for Fermi')
    fermi.add_argument('--fermi', type=str, dest='fermi_path',
        help='Path to Fermi wrapper', default='runfermi.pl')
    fermi.add_argument('--fermi_param', type=str, dest='fermi_param',
        help='Path to Fermi executable', default='fermi')                       
    bowtie = pbrazi.add_argument_group('Arguments for Bowtie2')
    bowtie.add_argument('--bowtie', type=str, dest='bowtie_path',
        help='Path to Bowtie2', default='bowtie2')
    samtools = pbrazi.add_argument_group('Arguments for Samtools')
    samtools.add_argument('--sam', type=str, dest='sam_path',
        help='Path to Samtools', default='samtools')
    bbduk = pbrazi.add_argument_group('Arguments for BBDuk')
    bbduk.add_argument('--bbduk', type=str, dest='bbduk_path',
        help='Path to bbduk executable', default='bbduk')
    bbduk.add_argument('--contaminant', type=str, dest='confile',
        help='Path to contaminant reference fasta file')
    bbduk.add_argument('--adapters', type=str, dest='adapters', nargs='+',
        help='Path to adapter and overrepresented sequence fasta file')
    opts = pbrazi.parse_args()

    prep = Prepper(opts.input_path)
    config = prep.prepInputs()
    print(config)
    main(opts.bbduk_path, opts.bowtie_path, opts.spades_path, opts.sga_path, opts.ngopt_path,
        opts.fermi_path, opts.abyss_path, opts.spades_param, opts.sga_param, opts.ngopt_param,
        opts.fermi_param, opts.abyss_param, config , opts.ref_path,
        opts.threads, opts.mode, opts.out_path, opts.adapters, opts.name)

