#!/usr/bin/env python
import argparse
import subprocess
import sys
sys.path.append('/net/shendure/vol1/home/cusanovi/bin/Levenshtein/')
import Levenshtein
import gzip
import io
import cStringIO
io_method = cStringIO.StringIO

parser = argparse.ArgumentParser(description='A program to fix erroneous barcodes in scATAC data.')
parser.add_argument('-F1','--Fastq1', help='Input fastq file read 1',dest='fastq1',required=True)
parser.add_argument('-F2','--Fastq2', help='Input fastq file read 2',dest='fastq2',required=True)
parser.add_argument('-H','--',help='Hash table',dest='Hash',required=True)
parser.add_argument('-O','--output', help='Output file ',dest='output',required=True)
#parser.add_argument('-L','--log', help='Output log file',dest='logfile',required=True)
args = parser.parse_args()

def submitter(commander):
    """Submits commands directly to the command line and waits for the process to finish."""
    submiting = subprocess.Popen(commander,shell=True)
    submiting.wait()

print "building a dictionaray of hashes used in experiment\n" 
hashList = open(args.Hash, "r")
hashDict = {}
for h in hashList.readlines():
    fields = h.strip().split()
    well, seq = fields[0], fields[1]
    hashDict[seq] = well

print "searching for hashIDs in fastq reads\n"
if1 = gzip.open(args.fastq1,'rb')
infastq1 = io.BufferedReader(if1)
if2 = gzip.open(args.fastq2,'rb')
infastq2 = io.BufferedReader(if2)
outF = open(args.output, "w")
for line in infastq2:
    seqbarc = line.strip('@').split(':')[0]
    read2 = infastq2.next().strip()
    plus2 = infastq2.next()
    qual2 = infastq2.next().strip()
    dumpbarc = infastq1.readline()
    read1 = infastq1.readline().strip()
    plus1 = infastq1.readline()
    qual1 = infastq1.readline().strip()

    try:
        hashDict[read2[0:10]]
        outF.write("%s\t%s\t%s\t%s\t%s\n" %("scichem3L",  seqbarc, read1[0:8], hashDict[read2[0:10]], 1))
    except KeyError:
        pass

outF.close()

