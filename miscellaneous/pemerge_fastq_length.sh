#$ -S /bin/bash
#$ -N pemerge_fastq_length
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/pemerge_length_fastq_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/pemerge_length_fastq_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=40G
#$ -q long.q@cbsubrito
#$ -t 1-8

WRK=/workdir/users/acv46/stool_PROSeq2/fastq/clean
export OMP_NUM_THREADS=8
export PATH=/programs/bwa-0.7.17:$PATH

LIST=$WRK/samples.txt
  DESIGN=$(sed -n "${SGE_TASK_ID}p" $LIST)
  NAME=`basename "$DESIGN"`

cd $WRK

echo "${NAME} -- merging PE reads"

bwa pemerge -m -t 8 \
	${NAME}_R1.fastq.gz ${NAME}_R2.fastq.gz \
	> ${NAME}_merged.fastq

bwa pemerge -u -t 8 \
        ${NAME}_R1.fastq.gz ${NAME}_R2.fastq.gz \
        > ${NAME}_unmerged.fastq

echo "${NAME} -- getting length distributions for fastqs"

awk '{if(NR%4==2) print length($1)}' ${NAME}_merged.fastq | \
	sort -n | uniq -c | \
	sed 's/^ *//g' | tr ' ' '\t' \
	> ${NAME}_merged_readlengths.txt

awk '{if(NR%4==2) print length($1)}' ${NAME}_unmerged.fastq | \
        sort -n | uniq -c | \
        sed 's/^ *//g' | tr ' ' '\t' \
        > ${NAME}_unmerged_readlengths.txt

echo "${NAME} -- gzipping fastqs"

gzip ${NAME}_merged.fastq & gzip ${NAME}_unmerged.fastq


