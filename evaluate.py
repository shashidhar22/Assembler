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
from collections import namedtuple
from operator import attrgetter
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

def ScafoldMetric(assemblies, assemblers, gsize):
    outfile = open('AssemblyMetrics.csv','w')
    outfile.write('Assembler\tScaffoldLength\tNG(%)\n')
    for assembly, assembler in zip(assemblies, assemblers):
        tlen = 0
        for contigs in FastaParser(assembly):
            tlen += contigs.length
            ngval = (tlen * 100)/float(gsize)
            outfile.write('{0}\t{1}\t{2:.2f}\n'.format(assembler, contigs.length, ngval))
    outfile.close()
    return
def ScafoldCounter(assemblies, assemblers, bsize):
    metric_matrix = list()
    metric_colnames = assemblers
    metric_rownames = np.linspace(0, 230000, bsize)
    contig_dict = dict()
    for assembly, assembler in zip(assemblies, assemblers):
        contig_list = list()
        for contigs in FastaParser(assembly):
            contig_list.append(contigs.length)
        contig_dict[assembler] = contig_list  #list(np.histogram(contig_list, metric_rownames)[0]) + [0]
    contig_table = pd.DataFrame(contig_dict, index=metric_rownames)
    return(contig_table)

if __name__ == '__main__':
    troch =  argparse.ArgumentParser(prog='Hummingbird')
    troch.add_argument('-f', '--assemblies', type=str, dest='assemblies', nargs='+',
        help='Assembly file list')
    troch.add_argument('-a', '--assemblers', type=str, dest='assemblers', nargs='+',
        help='Assembler name list')
    troch.add_argument('-b', '--bsize', type=int, dest='bsize', help='Number of bins to create')
    troch.add_argument('-g', '--gsize', type=int, dest='bsize', help='Estimated genome size')
    opts = troch.parse_args()
    #contig_table = ScafoldCounter(opts.assemblies, opts.assemblers, opts.bsize)
    #contig_table.to_csv('AssemblyStats.csv',header=True, sep=',')
    ScafoldMetric(opts.assemblies, opts.assemblers, opts.gsize)
