#$ -S /bin/bash
#$ -N filter_mapq
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/filter_mapq_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/filter_mapq_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=20G
#$ -q long.q@cbsubrito
#$ -t 1-16

# filter bam by mapq >= 20
# assumes each bam file will have at least one read with each possible mapq score (0-60)

WRK=/workdir/users/acv46/stool_PROSeq2/deseq/bam

export OMP_NUM_THREADS=8
export PATH=/programs/samtools-1.11/bin:$PATH

LIST=$WRK/bams.txt
  DESIGN=$(sed -n "${SGE_TASK_ID}p" $LIST)
  NAME=`basename "$DESIGN"`

# initialize text file
# printf '%s\n' {0..60} | sed "1 i\MAPQ" > mapq_scores.txt

BAM=$WRK/${NAME}.bam

samtools view -@ ${OMP_NUM_THREADS} $BAM | \
	awk -F$'\t' '$5 >= 20' | \
	> $WRK/${NAME}_MAPQ20.sam

lines=$(samtools  | wc -l)

if [[ "$lines" == 62 ]]; then

	echo "${NAME}_mapq.tmp has ${lines} lines -- GOOD!"
	paste $WRK/mapq_scores.txt $WRK/${NAME}_mapq.tmp \
	> tmpout && mv tmpout $WRK/mapq_scores.txt

else

	echo "${NAME}_mapq.tmp has ${lines} lines -- BAD! Must be 62. Skipping..."

fi

rm $WRK/${NAME}_mapq.tmp
