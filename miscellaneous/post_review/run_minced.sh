#$ -S /bin/bash
#$ -N minced
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq3/log/minced_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq3/log/minced_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq3
#$ -l h_vmem=50G
#$ -q long.q@cbsubrito
#$ -t 1-2
#$ -tc 2
#$ -pe parenv 4

# runs minced
# installed via makefile

# run on subset of all contigs >= 1000 nt

WRK=/workdir/users/acv46/stool_PROSeq3
#BLASTDB=$WRK/phage/virsorter
MIN=/home/acv46/minced/minced
CONTIG=$WRK/deseq/ref

export OMP_NUM_THREADS=4

DESIGN_FILE=$WRK/assembly/names.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`echo ${DESIGN} | sed "s/_25May2022//g"`

CONTIG=$WRK/transcriptomes/${NAME}/index/contigs_1000.fasta

OUT=$WRK/annotation/${NAME}_minced
mkdir -p $OUT
cd $OUT

$MIN -gffFull -spacers \
	-minNR 3\
	-minRL 20 \
	-maxRL 50 \
	-minSL 22 \
	-maxSL 55 \
	$CONTIG \
	${NAME}.txt \
	${NAME}.gff

