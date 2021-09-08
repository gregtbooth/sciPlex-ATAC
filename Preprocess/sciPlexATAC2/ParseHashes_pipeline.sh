# Below are the commands that were used to parse raw hash sequencing data into a count tables.
# The output is to be used to assess whether hash swapping between cells has occurred 
# under various sample prep conditions.
#-----------------------
# demultiplex my hash sequences which were spiked into sequencing run
#-----------------------

# run deMultiplexHashes.sh script to demultiplex sequencing output to fastq files.  

#----------------------
# Run following commands in directory containing the demultiplexed fastq files
#----------------------

WORKING_DIR=`pwd`
cd $WORKING_DIR

#
# Adjust these paths as necessary


SCRIPTS_DIR=~/bin
DATAMASH_PATH=~/datamash # Datamash is required and can be downloaded here: https://www.gnu.org/software/datamash/
HASH_BARCODES_FILE=~/bin/SpacePlate19and1.txt

# RT_BARCODES_FILE is a two column tab-delimited file where
# column 1 = well id
# column 2 = RT barcode

RT_BARCODES_FILE=~/bin/RT_indices_all

# If you have only one biological sample,
# just write the sample name to a file called combinatorial.indexing.key

SAMPLE_NAME="scichem2"
echo "$SAMPLE_NAME" >$WORKING_DIR/combinatorial.indexing.key

# Note that I run these processes using our cluster, which requires resource requests and specification of modules to load.  

# In various steps of the pipeline,
# one qsub job will be submitted for each BATCH_SIZE PCR wells of the experiment.
#
# If you have >1 PCR plate, you may want to increase BATCH_SIZE
# so you don't take up too many resources.
# The most memory intensive script that runs in batches in this manner
# uses 10 GB per qsub, or 240 GB total for 96 PCR wells and BATCH_SIZE=4
#

BATCH_SIZE=10

module unload mpc/0.8.2
module unload mpfr/3.1.0
module unload gmp/5.0.2
module load gcc/8.1.0
module load R/3.5.1 

#-------------------------------------------------------------------------------
# Put read 1 info (RT well, UMI) into read 2 read name
#-------------------------------------------------------------------------------

cd $WORKING_DIR

mkdir combined-fastq
mkdir file-lists-for-r1-info-munging
mkdir put-r1-info-in-r2-logs

ls fastq/ | grep _R1_ | grep -v Undetermined | split -l $BATCH_SIZE -d - file-lists-for-r1-info-munging/

ls file-lists-for-r1-info-munging | while read BATCH; do
    qsub -P trapnelllab $SCRIPTS_DIR/put-read1-info-in-read2.sh \
        $WORKING_DIR/fastq                                      \
        $WORKING_DIR/file-lists-for-r1-info-munging/$BATCH      \
        $SCRIPTS_DIR/                                           \
        $RT_BARCODES_FILE                                       \
        $WORKING_DIR/combinatorial.indexing.key                 \
        $WORKING_DIR/combined-fastq                             \
        $WORKING_DIR/put-r1-info-in-r2-logs
done

#-------------------------------------------------------------------------------
# Parse Hash Barcodes
#-------------------------------------------------------------------------------

cd $WORKING_DIR

mkdir file-lists-for-trimming
mkdir hashed-fastq
mkdir hashed-logs

ls combined-fastq/ | split -l $BATCH_SIZE -d - file-lists-for-trimming/

ls file-lists-for-r1-info-munging | while read BATCH; do
    qsub -P trapnelllab $SCRIPTS_DIR/parse_hash.sh \
        $WORKING_DIR/combined-fastq                             \
        $WORKING_DIR/file-lists-for-trimming/$BATCH             \
        $SCRIPTS_DIR/                                    \
        $HASH_BARCODES_FILE                                     \
        $WORKING_DIR/combinatorial.indexing.key                 \
        $WORKING_DIR/hashed-fastq                               \
        $WORKING_DIR/hashed-logs
done


#-------------------------------------------------------------------------------
# Generate Hash Tables 
#-------------------------------------------------------------------------------


cd $WORKING_DIR
mkdir $WORKING_DIR/hashRDS

zcat hashed-fastq/*.gz \
    | $DATAMASH_PATH -g 2,3,4 count 3 \
    | $DATAMASH_PATH -g 1 sum 4 \
    | awk -v S=$SAMPLE_NAME '{OFS="\t";} {print S, $0}' \
    > $WORKING_DIR/hashReads.per.cell


zcat hashed-fastq/*.gz \
    | uniq \
    | $DATAMASH_PATH -g 2,3,4 count 3 \
    | $DATAMASH_PATH -g 1 sum 4 \
    | awk -v S=$SAMPLE_NAME '{OFS="\t";} {print S, $0}' \
    > $WORKING_DIR/hashUMIs.per.cell


Rscript $SCRIPTS_DIR/knee-plot.R            \
    $WORKING_DIR/hashUMIs.per.cell          \
    $WORKING_DIR/hashed-logs


zcat hashed-fastq/*.gz \
    | uniq \
    | $DATAMASH_PATH -g 1,2,4,5 count 3  \
    > $WORKING_DIR/hashRDS/hashTable.out 

paste $WORKING_DIR/hashUMIs.per.cell  $WORKING_DIR/hashReads.per.cell \
    | cut -f 1,2,6,3 \
    | awk 'BEGIN {OFS="\t";} {dup = 1-($3/$4); print $1,$2,$3,$4,dup;}' \
    > $WORKING_DIR/hashDupRate.txt







