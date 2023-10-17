#$ -S /bin/bash
#$ -N bam2cov
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq3/log/bam2cov_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq3/log/bam2cov_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq3
#$ -l h_vmem=30G
#$ -q long.q@cbsubrito
#$ -t 1-22
#$ -tc 4
#$ -pe parenv 4

# split presorted bams into forward and reverse
# generate indices
# generate coverage reports with bedtools

export OMP_NUM_THREADS=4
export PATH=/programs/samtools-1.15.1-r/bin:$PATH
export PATH=/programs/bedtools2-2.29.2/bin:$PATH
WRK=/workdir/users/acv46/stool_PROSeq3/transcriptomes

List=${WRK}/bams.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" ${List})
        NAME=`basename "$DESIGN"`

samp=$(echo "${NAME}" | awk -F_ '{print $1}')
base=$(echo "${NAME}" | sed "s/_dedup_QC_end.sort.bam//g")
type=$(echo "${NAME}" | awk -F_ '{print $3}')

BAM=${WRK}/${samp}/old_bam/${NAME}

mkdir -p ${WRK}/${samp}/coverage
cd ${WRK}/${samp}/coverage

if [ "${type}" == "PROseq" ]; then

  echo "${NAME} -- starting pipeline for PROseq library"

  samtools view -@ ${OMP_NUM_THREADS} -f 147 -F 3852 -q 30 -b ${BAM} > ${base}_rev1.bam
  samtools view -@ ${OMP_NUM_THREADS} -f  99 -F 3852 -q 30 -b ${BAM} > ${base}_rev2.bam
  samtools view -@ ${OMP_NUM_THREADS} -f 131 -F 3868 -q 30 -b ${BAM} > ${base}_fwd1.bam
  samtools view -@ ${OMP_NUM_THREADS} -f  67 -F 3884 -q 30 -b ${BAM} > ${base}_fwd2.bam

  samtools merge \
    -@ ${OMP_NUM_THREADS} \
    -o ${base}_fwd_tosort.bam \
    ${base}_fwd1.bam ${base}_fwd2.bam
  samtools sort \
    --write-index \
    -@ ${OMP_NUM_THREADS} \
    -o ${base}_fwd.bam \
    ${base}_fwd_tosort.bam

  samtools merge \
    -@ ${OMP_NUM_THREADS} \
    -o ${base}_rev_tosort.bam \
    ${base}_rev1.bam ${base}_rev2.bam
  samtools sort \
    --write-index \
    -@ ${OMP_NUM_THREADS} \
    -o ${base}_rev.bam \
    ${base}_rev_tosort.bam

  rm ${base}_fwd_tosort.bam ${base}_rev_tosort.bam \
     ${base}_fwd1.bam ${base}_fwd2.bam \
     ${base}_rev1.bam ${base}_rev2.bam

  bedtools genomecov -d -ibam ${base}_fwd.bam > ${base}_fwd_depth.txt
  bedtools genomecov -d -ibam ${base}_rev.bam > ${base}_rev_depth.txt

  echo "${NAME} -- finished pipeline for PROseq library"

elif [ "${type}" == "TEXp" ] || \
     [ "${type}" == "TEXm" ]; then

  echo "${NAME} -- starting pipeline for TEX library"

  samtools view -@ ${OMP_NUM_THREADS} -f 147 -F 3852 -q 30 -b ${BAM} > ${base}_fwd1.bam
  samtools view -@ ${OMP_NUM_THREADS} -f  99 -F 3852 -q 30 -b ${BAM} > ${base}_fwd2.bam
  samtools view -@ ${OMP_NUM_THREADS} -f 131 -F 3868 -q 30 -b ${BAM} > ${base}_rev1.bam
  samtools view -@ ${OMP_NUM_THREADS} -f  67 -F 3884 -q 30 -b ${BAM} > ${base}_rev2.bam

  samtools merge \
    -@ ${OMP_NUM_THREADS} \
    -o ${base}_fwd_tosort.bam \
    ${base}_fwd1.bam ${base}_fwd2.bam
  samtools sort \
    --write-index \
    -@ ${OMP_NUM_THREADS} \
    -o ${base}_fwd.bam \
    ${base}_fwd_tosort.bam

  samtools merge \
    -@ ${OMP_NUM_THREADS} \
    -o ${base}_rev_tosort.bam \
    ${base}_rev1.bam ${base}_rev2.bam
  samtools sort \
    --write-index \
    -@ ${OMP_NUM_THREADS} \
    -o ${base}_rev.bam \
    ${base}_rev_tosort.bam

  rm ${base}_fwd_tosort.bam ${base}_rev_tosort.bam \
     ${base}_fwd1.bam ${base}_fwd2.bam \
     ${base}_rev1.bam ${base}_rev2.bam

  bedtools genomecov -d -ibam ${base}_fwd.bam > ${base}_fwd_depth.txt
  bedtools genomecov -d -ibam ${base}_rev.bam > ${base}_rev_depth.txt

  echo "${NAME} -- finished pipeline for TEX library"

else

  echo "${NAME} -- error parsing library type"

fi
