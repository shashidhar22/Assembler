import os
import sys
import csv
import time
import logging
from collections import namedtuple


class Fastq:

    def __init__(self, fastqfile, outdir, phred):
        self.fastq = fastqfile
        self.outdir = outdir
        self.phred = phred
        FORMAT = '%(asctime)-15s : %(levelname)-8s :  %(message)s'
        logging.basicConfig(format=FORMAT)
        self.phreddict = self.phredmap()
        

    def phredmap(self):
        phreddict = dict()
        if self.phred == "phred33":
            for asciis, quals in zip(range(33,126), range(0,92)):
                phreddict[asciis] = quals
        return(phreddict)

    def formatChecker(self, header, seq, sheader, quals, line, phred):
        if header[0] != "@":
            return(1)
        if sheader == "" or sheader[0] != "+":
            return(2)
        if len(seq) != len(quals):
            return(3)
        if seq == '':
            return(4)

    def qualmasker(self, seq, quals):
        maskedseq = list()
        for base, qual in zip(seq, quals):
            if qual < 20:
                maskedseq.append('N')
            else:
                maskedseq.append(base)
        return(maskedseq, quals)
        
    def read(self):
        fqhandle = open(self.fastq, 'r')
        Record = namedtuple('Record',['header', 'sheader', 'seq', 'quals'])
        lineno = 0
        while True:
            try:
                #Get fastq record
                header = next(fqhandle).strip()
                lineno
                seq = [base for base in next(fqhandle).strip()]
                sheader = next(fqhandle).strip()
                quals = [self.phreddict[int(ord(qual))] for qual in next(fqhandle).strip()]
                lineno += 1
            
                #Check if record adheres to fastq format
                check = self.formatChecker(header, seq, sheader, quals, lineno, self.phred)
                if check == 1:
                    logging.error('Invalid header in fastq read ; record number : {1}'.format(lineno))
                    raise NameError('Invalid header in fastq read; record number : {0}'.format(lineno))
                if check == 2:
                    logging.error('Invalid secondary header in fastq read ; record number : {1}'.format(lineno))
                    raise NameError('Invalid secondary header in fastq read; record number : {0}'.format(lineno))
                if check == 3:
                    logging.error('Sequence and quality strings of unequl length in fastq read; record number : {0}'.format(lineno))
                    raise NameError('Sequence and quality strings of unequal length in fastq read; record number : {0}'.format(lineno))
                if check == 4:
                    logging.error('Sequence data missing; record number : {0}'.format(lineno))
                    raise NameError('Sequence data missing; record number : {0}'.format(lineno))

                #Optionally mask low quality reads
                #seq, quals = self.qualmasker(seq, quals)

                #Return record
                record = Record(header, sheader, seq, quals)
                yield(record)
            except StopIteration:
                logging.info('End of file')
                break
            except NameError:
                break

        return



class Pileup:
    #Class to parse the pileup format and extract relevant information
    
    def __init__(self, pileupfile, outdir):
        self.pileup = pileupfile
        self.outdir = outdir
        FORMAT = '%(asctime)-15s : %(levelname)-8s :  %(message)s'
        logging.basicConfig(format=FORMAT)

    def getAlt(self, alt, qual):
        alt = list()
        altdir = list()
        altevi = dict()
        altqual = list()
        for base in alt:
            if base == '.' or base == ',':
                al
                
        
    def read(self):
        phandle = open(self.pileup)
        Baseinfo = namedtuple("Baseinfo", ['scaf', 'pos', 'ref', 'cov', 'alt', 'altdir', 'altevi', 'altqual'])
        while True:
            try:
                lines = lines.split('\t').strip()
                scaf = lines[0]
                pos = lines[1]
                ref = lines[2]
                cov = lines[3]
                alt = lines[4]
                qual = lines[5]
            
