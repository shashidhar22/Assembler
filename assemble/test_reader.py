#! /usr/local/bin/python3.4
from readers import Fastq
test = Fastq('../fq/test_reader.fastq','./','phred33')
reader = test.read()
for records in reader:
    print(records) 
