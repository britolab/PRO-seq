#$ -S /bin/bash
#$ -N bmtagger
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq2/log/bmtagger_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq2/log/bmtagger_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=30G
#$ -q long.q@cbsubrito
#$ -t 1-8

#Set directories
WRK=/workdir/users/acv46/stool_PROSeq2
IN=$WRK/fastq/raw

#Create design file of file names
DESIGN_FILE=$WRK/scripts/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

## This script runs BMTagger which removes human reads from metagenomic data

################# BMTAGGER ######################

export PATH=/programs/ncbi-blast-2.3.0+/bin:$PATH
BMTAGGER=/programs/bmtools/bmtagger

REFGENOME=/home/britolab/refdbs/HUMAN/Homo_sapiens_assembly19.fasta
CONFIG=/workdir/users/acv46/scripts/cleanReads/bmtagger.conf

gunzip $IN/${NAME}_R1.fastq.gz
gunzip $IN/${NAME}_R2.fastq.gz

echo ${NAME} bmtagger start

${BMTAGGER}/bmtagger.sh \
-C $CONFIG \
-b ${REFGENOME}.bitmask \
-x ${REFGENOME}.srprism \
-T $IN -q1 \
-1 $IN/${NAME}_R1.fastq \
-2 $IN/${NAME}_R2.fastq \
-o $WRK/fastq/clean/${NAME} \
-X

echo $NAME bmtagger complete

gzip $IN/${NAME}_R1.fastq
gzip $IN/${NAME}_R2.fastq
