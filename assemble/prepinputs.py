import os
import sys
import re
import glob
from collections import namedtuple
from assemble.readers import Fastq


class Identifier:

    def __init__(self, record):
        self.rec = record


    def isIlluminaOld(self):
        header_regex = re.compile('@\w+-?\w+:\d+:\d+:\d+:\d+#\d*')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)


    def isIlluminaNew(self):
        #@D00468:24:H8ELMADXX:1:1101:1470:2237 1:N:0:2
        header_regex = re.compile('@\w+-?\w+:\d+:\w+-?\w+:\d+:\d+:\d+:\d+\s\d:\w+:\w+:\w*')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)
  

    def isSraOld(self):
        #@SRR037455.1 HWI-E4_6_30ACL:4:1:0:29 length=35
        #@SRR902931.1 HWI-ST1384:61:D1DJ4ACXX:8:1101:1240:2015 length=50
        header_regex = re.compile('@\w+\.?\w? \w+-\w+:\d+:\d+:\d+:\d+ length=\d+')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)

    def isSraNew(self):
        header_regex = re.compile('@\w+\.?\w? \w+-\w+:\d+:\w+:\d+:\d+:\d+:\d+ length=\d+')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)

    def isPacbio(self):
        #@m160113_152755_42135_c100906712550000001823199104291667_s1_p0/15/7044_26271
        header_regex = re.compile('@\w+/\d+/\d+_\d+')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)

class Metrics:

    def __init__(self, fastq):
        self.fastq = fastq

    def avgReadLen(self):
        fastq_reader = Fastq(self.fastq, './', 'phred33')
        total_length = 0
        total_reads = 0
        for lines in fastq_reader.read():
            total_length += len(lines.seq)
            total_reads += 1
            if total_reads >= 100000:
                break

        avg_length = total_length/float(total_reads)
        return(avg_length)

def main(path):
    path = os.path.abspath(path)
    print(path)
    files = glob.glob('{0}/*.fastq*'.format(path))
    experiment = dict()
    for fastq in files:
        reader = Fastq(fastq, './', 'phred33')
        Sample = namedtuple('Sample', ['sample', 'library', 'files', 'prep', 'paired'])
        rec = next(reader.read())
        identifier = Identifier(rec)    
        metric = Metrics(fastq)
        isIllOld =  identifier.isIlluminaOld()
        isIllNew =  identifier.isIlluminaNew()
        isSraOld = identifier.isSraOld()
        isSraNew = identifier.isSraNew()
        isPac = identifier.isPacbio()
        seqType = ''
        libType = ''
        if isIllOld:
            paired_regex = re.compile('@\w+-?\w+:\d+:\d+:\d+:\d+#\d')
            lib = re.findall(paired_regex, rec.header)[0]
            paired = False
            seqType = 'Illumina'
            if metric.avgReadLen():
                libType = 'Short'
        elif isIllNew:
            paired_regex = re.compile('@\w+-?\w+:\d+:\w+-?\w+:\d+:\d+:\d+:\d+\s')
            lib = re.findall(paired_regex, rec.header)[0]
            paired = False
            seqType = 'Illumina'
            if metric.avgReadLen():
                libType = 'Short'
        elif isSraOld:
            paired_regex = re.compile('@\w+\.?\w? \w+-\w+:\d+:\w+:\d+:\d+:\d+:\d+ length=\d+')
            lib = re.findall(paired_regex, rec.header)[0]
            paired = False
            seqType = 'Illumina'
            if metric.avgReadLen():
                libType = 'Short'
        elif isSraNew:
            paired_regex = re.compile('@\w+\.?\w? \w+-\w+:\d+:\w+:\d+:\d+:\d+:\d+ length=\d+')
            lib = re.findall(paired_regex, rec.header)[0]
            paired = False
            seqType = 'Illumina'
            if metric.avgReadLen():
                libType = 'Short'
        elif isPac:
            lib_regex = re.compile('@\w+_\d+_\d+_\w+')
            lib = re.findall(lib_regex, rec.header)[0]
            paired = False
            seqType = 'Pacbio'
            if metric.avgReadLen():
                libType = 'Long'
        try:
            paired = True
            experiment[lib] = Sample(lib, seqType, [experiment[lib].files[0], fastq], libType, paired)
        except KeyError:            
            experiment[lib] = Sample(lib, seqType, [fastq], libType, paired)

    print(experiment)
    return(experiment)

if __name__ == '__main__':
    path = os.path.abspath(sys.argv[1])
    main(path)

    
