#$ -S /bin/bash
#$ -N merge_filter_pileup
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq3/log/mfp_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq3/log/mfp_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq3
#$ -l h_vmem=120G
#$ -q long.q@cbsubrito
#$ -t 1-4
#$ -tc 2
#$ -pe parenv 8

# merge replicate bam files
# filter merged file
# create bedtools pileup

export OMP_NUM_THREADS=8
export PATH=/programs/bedtools2-2.29.2/bin:$PATH
export PATH=/programs/samtools-1.15.1-r/bin:$PATH
WRK=/workdir/users/acv46/stool_PROSeq3/transcriptomes
OUT=$WRK/pileup
mkdir -p $OUT

qscore=30

List=$WRK/merged_US2.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $List)
        NAME=`basename "$DESIGN"`

echo "${NAME} -- starting MFP pipeline"
samp=$(echo ${NAME} | awk -F _ '{print $1}')
type=$(echo ${NAME} | awk -F _ '{print $2}')

######################
## MERGE and FILTER ##
######################

# merge and resort by coordinates
# filter for alignments with following characteristics
## mapq score >= 30
## read paired and in proper pair
## read not unmapped, mate not unmapped, not not primary alignment, not supplementary alignment
# -f 3 -F 2316

cd ${OUT}
mgbam=${OUT}/${NAME}_merged_filtered_q${qscore}.bam

if [ ! -f ${mgbam} ]; then

  samtools merge \
    -@ $OMP_NUM_THREADS \
    -n - \
    $WRK/${samp}/${samp}_*_${type}_dedup_QC_end.sort.bam | \
  samtools sort \
    -@ $OMP_NUM_THREADS | \
  samtools view \
    -@ $OMP_NUM_THREADS \
    -f 3 -F 2316 \
    -b -q ${qscore} \
    -o ${mgbam}

fi

echo "${NAME} -- bam files merged and filtered"

############
## PILEUP ##
############

# see /workdir/users/acv46/EC_PROSeq/pileup/_README
# P5 Illumina adapter on 3' end of nascent transcript
# 5' end of first in pair gives 3' end of read
# 5' end of second in pair gives 5' end of read
# flip strands

# sort the files before combining, dumbass

if [ ${type} = "PROseq" ]; then

  echo "${NAME} -- generating stranded paired-end coverage for PRO-seq alignments"

  bedtools genomecov -pc -d -strand + -ibam ${mgbam} | \
    sort -nk1.6 -nk2 \
    > ${NAME}_minus_q${qscore}_full.cov

  echo " --> ${NAME}_minus_q${qscore}_full.cov -- $(wc -l < ${NAME}_minus_q${qscore}_full.cov)"

  bedtools genomecov -pc -d -strand - -ibam ${mgbam} | \
    sort -nk1.6 -nk2 \
    > ${NAME}_plus_q${qscore}_full.cov

  echo " --> ${NAME}_plus_q${qscore}_full.cov -- $(wc -l < ${NAME}_plus_q${qscore}_full.cov)"

  samtools view -@ $OMP_NUM_THREADS -f 64 -b ${mgbam} | \
    bedtools genomecov -5 -d -strand + -ibam - | \
    sort -nk1.6 -nk2 \
    > ${NAME}_minus_q${qscore}_3p.cov

  echo " --> ${NAME}_minus_q${qscore}_3p.cov -- $(wc -l < ${NAME}_minus_q${qscore}_3p.cov)"

  samtools view -@ $OMP_NUM_THREADS -f 64 -b ${mgbam} | \
    bedtools genomecov -5 -d -strand - -ibam - | \
    sort -nk1.6 -nk2 \
    > ${NAME}_plus_q${qscore}_3p.cov

  echo "${NAME}_plus_q${qscore}_3p.cov -- $(wc -l < ${NAME}_plus_q${qscore}_3p.cov)"

  # because proper pairs map to opposite strands, 5' fragment end is correct strand

  samtools view -@ $OMP_NUM_THREADS -f 128 -b ${mgbam} | \
    bedtools genomecov -5 -d -strand + -ibam - | \
    sort -nk1.6 -nk2 \
    > ${NAME}_plus_q${qscore}_5p.cov

  echo "${NAME}_plus_q${qscore}_5p.cov -- $(wc -l < ${NAME}_plus_q${qscore}_5p.cov)"

  samtools view -@ $OMP_NUM_THREADS -f 128 -b ${mgbam} | \
    bedtools genomecov -5 -d -strand - -ibam - | \
    sort -nk1.6 -nk2 \
    > ${NAME}_minus_q${qscore}_5p.cov

  echo "${NAME}_minus_q${qscore}_5p.cov -- $(wc -l < ${NAME}_minus_q${qscore}_5p.cov)"

elif [ ${type} = "TEXm" ] || [ ${type} = "TEXp" ] || [ ${type} = "RNAseq" ]; then

  echo "${NAME} -- generating stranded paired-end coverage for RNA-seq alignments"

  bedtools genomecov -pc -d -strand + -ibam ${mgbam} | \
    sort -nk1.6 -nk2 \
    > ${NAME}_plus_q${qscore}_full.cov

  echo "${NAME}_plus_q${qscore}_full.cov -- $(wc -l < ${NAME}_plus_q${qscore}_full.cov)"

  bedtools genomecov -pc -d -strand - -ibam ${mgbam} | \
    sort -nk1.6 -nk2 \
    > ${NAME}_minus_q${qscore}_full.cov

  echo "${NAME}_minus_q${qscore}_full.cov -- $(wc -l < ${NAME}_minus_q${qscore}_full.cov)"

fi

#paste ${NAME}_*_q${qscore}_*.cov | \
#	awk '{print $1,$2,$3,$6,$9,$12,$15,$18,$21,$24}' | \
#	tr [[:blank:]] "\t" \
#	> ${NAME}_coverage_q${qscore}_correct3.txt

echo "${NAME} -- done!"

# remove intermediate files
#rm ${NAME}_*_q${qscore}_*.cov
#rm ${NAME}_PRO_merge_q${qscore}.bam
#rm ${NAME}_RNA_merge_q${qscore}.bam

