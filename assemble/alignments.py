import os
import sys
import pysam



def readFasta(fasta):
    fasta_handle = open(fasta)
    fasta_content = fasta_handle.read().split('>')[1:]
    fastaid = 1
    Fasta = namedtuple('Fasta',['header','seq','fid','length'])
    for lines in fasta_content:
        header = lines.split('\n')[0].split(' ')[0]
        sequence = ''.join(lines.split('\n')[1:]).upper()
        length = len(sequence)
        record = Fasta(header, sequence, fastaid, length)
        fastaid += 1
        yield record

def bamMetric(assembly, alignment):
    samfile = pysam.AlignmentFile(alignment,'rb')
    fasta = readFasta(assembly)
    
