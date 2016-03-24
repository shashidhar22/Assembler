import os
import sys
import Bio
import glob
import argparse
import itertools
import subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO
from itertools import repeat
from multiprocessing import Pool

def FastaParser(fasta):
    '''Parse fasta and return iterator with header and sequence'''
    fasta_handle = open(fasta)
    fasta_content = fasta_handle.read().split('>')[1:]
    reference = list()
    transcriptid = 0
    Fasta = namedtuple('Fasta',['header','seq','tid','length'])
    for lines in fasta_content:
        header = lines.split('\n')[0]
        sequence = ''.join(lines.split('\n')[1:]).upper()
        length = len(sequence)
        record = Fasta(header, sequence, transcriptid, length)
        transcriptid += 1
        yield record

def ScafoldCounter(assemblies, assemblers):
    metric_matrix = list()
    metric_colnames = assemblers
    metric_rownames = np.linspace(0, 230000, 1000)
    contig_dict = dict()
    for assembly, assembler in zip(assemblies, assemblers):
        contig_list = list()
        for contigs in FastaParser(assembly):
            contig_list.append(contigs.length)
        contig_dict[assembler] = np.histogram(contig_list, metric_rownames)
    contig_table = pd.DataFrame(contig_dict, index=metric_rownames)
    print(contig_table)
    return(contig_table)

if __name__ == '__main__':
    troch =  argparse.ArgumentParser(prog='Hummingbird')
    troch.add_argument('-f', '--assemblies', type=str, dest='assemblies', nargs='+',
        help='Assembly file list')
    troch.add_argument('-a', '--assemblers', type=str, dest='assemblers', nargs='+',
        help='Assembler name list')
    opts = troch.parse_args()
    contig_dict = ScafoldCounter(opts.assemblies, opts.assemblers)
