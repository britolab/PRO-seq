#$ -S /bin/bash
#$ -N minced
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq2/log/minced_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq2/log/minced_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=50G
#$ -q long.q@cbsubrito
#$ -t 1-2

# runs minced
# installed via makefile

# run on subset of all contigs >= 500 nt

WRK=/workdir/users/acv46/stool_PROSeq2
#BLASTDB=$WRK/phage/virsorter
MIN=/home/acv46/minced/minced
OUT=$WRK/phage/minced_out/standalone
CONTIG=$WRK/deseq/ref

export OMP_NUM_THREADS=8

DESIGN_FILE=$OUT/../../allcontigs2.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

OUT1=$OUT/${NAME}
mkdir -p $OUT1
cd $OUT1

$MIN -gffFull -spacers \
	-minNR 3\
	-minRL 20 \
	-maxRL 50 \
	-minSL 22 \
	-maxSL 55 \
	$CONTIG/${NAME}_allcontigs.short_500.fasta \
	${NAME}.txt \
	${NAME}.gff

