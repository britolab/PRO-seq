#$ -S /bin/bash
#$ -N contigcounts
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq2/log/contigcounts_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq2/log/contigcounts_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=40G
#$ -q long.q@cbsubrito
#$ -t 1-4

WRK=/workdir/users/acv46/stool_PROSeq2/deseq
BAM=$WRK/bam
FCT=$WRK/featureCounts/US2_goodbins
export OMP_NUM_THREADS=6

DESIGN_FILE=$FCT/US2_goodbins_bams.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

base=$(echo $NAME | awk -F. '{print $1}')
samtools view -F 260 $BAM/${NAME} | \
awk -F$'\t' '{print $3}' | \
sort | uniq -c > $FCT/${base}.primarymapped.txt

