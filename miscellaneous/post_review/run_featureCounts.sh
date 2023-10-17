#$ -S /bin/bash
#$ -N featureCounts
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq3/log/featureCounts_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq3/log/featureCounts_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq3
#$ -l h_vmem=100G
#$ -q long.q@cbsubrito
#$ -t 1-2
#$ -tc 2
#$ -pe parenv 12

# iterates over different GTF-bam pairings
# runs REV mapping (-s 2) for both RNA-seq and PRO-seq

WRK=/workdir/users/acv46/stool_PROSeq3/transcriptomes
FCT=/home/britolab/acv46/subread-2.0.2-source/bin/featureCounts
export OMP_NUM_THREADS=12

DESIGN_FILE=$WRK/base.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

annot=/workdir/users/acv46/stool_PROSeq3/annotation/${NAME}_bins_bakta/${NAME}_bins_merged.gff3

# for gff in ./*/*.gff3; do sed -n '/>/q;p' $gff >> ${samp}_bins_merged.gff3; done

OUT=$WRK/featureCounts/${NAME}
mkdir -p $OUT

# paired reads (-p)
# only count read pairs where both ends are aligned (-B)
# do not count read pairs if reads map to different contigs or different strands (-C)
# perform read counting at feature level (-f)

cd $OUT

echo "${NAME} -- filtering bams by MAPQ, etc."

# subset gff for contigs in bam

for bam in $WRK/${NAME}/*_dedup_QC_end.sort.bam; do
  grep "${bam}"
  bname=$(basename ${bam} | sed "s/_dedup_QC_end.sort.bam//g")
  samtools view \
    -@ $OMP_NUM_THREADS \
    -f 3 -F 2316 \
    -b -q 30 \
    -o $OUT/${bname}_q30.bam \
    ${bam}
  samtools view  \
    -@ $OMP_NUM_THREADS \
    $OUT/${bname}_q30.bam | \
    awk -F$"\t" '{print $3}' | sort | uniq \
    > ${bname}_bamcontigs.txt
done

cat ${OUT}/*_bamcontigs.txt | sort | uniq \
  > ${OUT}/${NAME}_contigs.txt

rm ${OUT}/*_bamcontigs.txt

grep -v "#" ${annot} | \
  awk -F$'\t' 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' \
  ${OUT}/${NAME}_contigs.txt - \
  > ${OUT}/${NAME}_bins.gff3

##samtools index $OUT/${NAME}_q30.bam $OUT/${NAME}_q30.bam.index

##samtools view -b $OUT/${NAME}_q30.bam \
##	${samp}_contigs.txt \
##	> $OUT/${NAME}_q30_sub.bam

# get sense (reverse) counts

echo "${NAME} -- getting sense feature counts"

$FCT \
  -a ${OUT}/${NAME}_bins.gff3 \
  -o $OUT/${NAME}_q30_fcounts_rev.txt \
  -F GTF \
  -t CDS,CRISPR,ncRNA,oriC,regulatory_region,rRNA,tmRNA,tRNA \
  -g ID \
  --extraAttributes Name,product \
  -s 2 \
  --countReadPairs \
  -f -p -B -C \
  -T $OMP_NUM_THREADS \
  --verbose \
  $OUT/*_q30.bam

# get antisense (forward) counts

echo "${NAME} -- getting antisense feature counts"

$FCT \
  -a ${OUT}/${NAME}_bins.gff3 \
  -o $OUT/${NAME}_q30_fcounts_fwd.txt \
  -F GTF \
  -t CDS,CRISPR,ncRNA,oriC,regulatory_region,rRNA,tmRNA,tRNA \
  -g ID \
  --extraAttributes Name,product \
  -s 1 \
  --countReadPairs \
  -f -p -B -C \
  -T $OMP_NUM_THREADS \
  --verbose \
  $OUT/*_q30.bam

echo "${NAME} -- removing intermediate bams"

rm $OUT/*_q30.bam
