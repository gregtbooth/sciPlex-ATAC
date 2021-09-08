#! usr/bin/bash

# script to request cluster time: 
qlogin -q trapnell-login.q -l mfree=15G -pe serial 1
# notes: mfree = memory needed,  serial = threads.  

# script for demultiplexing sequencing output. 
module load bcl2fastq/2.18

SEQPATH=~/NEXTSEQ/190518_NS500488_0826_AHJ5LCBGXB
SAMPLESHEET=/bin//HashSampleSheetF0102030405_nextera.csv


bcl2fastq -R $SEQPATH -o fastq/ --no-lane-splitting --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0 --sample-sheet $SAMPLESHEET

