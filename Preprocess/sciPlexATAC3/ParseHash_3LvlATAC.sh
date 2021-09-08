# Below are the commands that were used to parse raw hash sequencing data into a count tables.
# The output is to be used to assess whether hash swapping between cells has occured 
# under various sample prep conditions.


#----------------------
# Run commands in directory containing the demultiplexed fastq files
#----------------------

qlogin -q trapnell-login.q -l mfree=50G -pe serial 1

WORKING_DIR=`pwd`
cd $WORKING_DIR

SCRIPTS_DIR=~/bin
DATAMASH_PATH=~/datamash
HASH_BARCODES_FILE=/bin/Space_Plate10.txt


# If you have only one biological sample,
# just write the sample name to a file called combinatorial.indexing.key

SAMPLE_NAME="MLR"
echo "$SAMPLE_NAME" >$WORKING_DIR/combinatorial.indexing.key

# uses 10 GB per qsub, or 240 GB total for 96 PCR wells and BATCH_SIZE=4


module unload mpc/0.8.2
module unload mpfr/3.1.0
module unload gmp/5.0.2
module load gcc/8.1.0
module load R/3.5.2 
module load python/2.7.3

#-------------------------------------------------------------------------------
# Parse Hash Barcodes from trimmed, barcode corrected fastq (read2) 
#-------------------------------------------------------------------------------
cd $WORKING_DIR
mkdir $WORKING_DIR/parse_hash_out
python $SCRIPTS_DIR/ParseHash_3L_ATAC_GB.py -F1 $WORKING_DIR/fastqs/Undetermined_S0_R1_001.fastq.gz.out.fq.gz -F2 $WORKING_DIR/fastqs/Undetermined_S0_R2_001.fastq.gz.out.fq.gz -H $HASH_BARCODES_FILE -O $WORKING_DIR/parse_hash_out/Parsed_Hash_Table

#-------------------------------------------------------------------------------
# Parse Hash Barcodes
#-------------------------------------------------------------------------------

cd $WORKING_DIR
mkdir $WORKING_DIR/hashRDS

cat $WORKING_DIR/parse_hash_out/Parsed_Hash_Table \
    | $DATAMASH_PATH -s -g 2,3,4 count 3 \
    | $DATAMASH_PATH -g 1 sum 4 \
    | awk -v S=$SAMPLE_NAME '{OFS="\t";} {print S, $0}' \
    > $WORKING_DIR/hashRDS/hashReads.per.cell

cat $WORKING_DIR/parse_hash_out/Parsed_Hash_Table \
    | uniq \
    | $DATAMASH_PATH -s -g 2,3,4 count 3 \
    | $DATAMASH_PATH -g 1 sum 4 \
    | awk -v S=$SAMPLE_NAME '{OFS="\t";} {print S, $0}' \
    > $WORKING_DIR/hashRDS/hashUMIs.per.cell

Rscript $SCRIPTS_DIR/knee-plot.R            \
    $WORKING_DIR/hashRDS/hashUMIs.per.cell          \
    $WORKING_DIR/hashRDS

cat $WORKING_DIR/parse_hash_out/Parsed_Hash_Table \
    | uniq \
    | $DATAMASH_PATH -s -g 1,2,4,5 count 3  \
    > $WORKING_DIR/hashRDS/hashTable.out 

paste $WORKING_DIR/hashRDS/hashUMIs.per.cell $WORKING_DIR/hashRDS/hashReads.per.cell\
    | cut -f 1,2,6,3 \
    | awk 'BEGIN {OFS="\t";} {dup = 1-($3/$4); print $1,$2,$3,$4,dup;}' \
    > $WORKING_DIR/hashRDS/hashDupRate.txt








