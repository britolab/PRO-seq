#$ -S /bin/bash
#$ -N mapq_length
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/mapq_length_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/mapq_length_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2  #Working directory
#$ -l h_vmem=20G
#$ -q long.q@cbsubrito
#$ -t 1-4

WRK=/workdir/users/acv46/stool_PROSeq2/deseq/bam/pileup
export OMP_NUM_THREADS=8

LIST=$WRK/bams.txt
  DESIGN=$(sed -n "${SGE_TASK_ID}p" $LIST)
  NAME=`basename "$DESIGN"`

OUT=$WRK/stats
mkdir -p $OUT
cd $WRK

base=$(echo ${NAME} | sed "s/_merge.bam//g")

echo "${NAME} -- calculating mapq distribution"

samtools view -@ 8 ${NAME} | \
	awk -F$"\t" '{print $5}' | \
	sort | uniq -c | sed 's/^ *//g' | \
	tr ' ' '\t'  | sort -k2n \
	> $OUT/${base}_mapq.txt

echo "${NAME} -- calculating mapped read length distribution"

samtools view -@ 8 -f 3 -F 12 ${NAME} | \
	cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | \
	sort | uniq -c | sed 's/^ *//g' | \
        tr ' ' '\t'  | sort -k2n \
	> $OUT/${base}_mappedReadLength.txt
