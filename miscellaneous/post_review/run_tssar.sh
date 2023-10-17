#$ -S /bin/bash
#$ -N run_tssar
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq3/log/tssar_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq3/log/tssar_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq3
#$ -l h_vmem=50G
#$ -q long.q@cbsubrito
#$ -t 1-2
#$ -tc 2
#$ -pe parenv 8

# merge replicate bam files for dRNA-seq reads
# filter merged file
# run tssar (plus v. minus)

export OMP_NUM_THREADS=8
export PATH=/programs/bedtools2-2.29.2/bin:$PATH
export PATH=/programs/samtools-1.15.1/bin:$PATH
WRK=/workdir/users/acv46/stool_PROSeq3/transcriptomes

qscore=30
nreads=1000
length=100000

List=$WRK/base.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $List)
        NAME=`basename "$DESIGN"`

BAM=$WRK/${NAME}
OUT=$WRK/tssar/${NAME}_PRO
SEQ=${BAM}/index/contigs_1000.fasta
mkdir -p $OUT

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
mgsam_p=${OUT}/${NAME}_merged_filtered_q${qscore}_PLUS.sam
mgsam_m=${OUT}/${NAME}_merged_filtered_q${qscore}_MINUS.sam
CLIST=${OUT}/${NAME}_contigs_${length}nt_${nreads}reads_contigs.txt

if [[ ! -f ${mgsam_p} ]] || \
   [[ ! -f ${mgsam_m} ]]; then

  echo "${NAME} -- merging and filtering bams"

  # PROseq as plus
  samtools merge \
    -@ $OMP_NUM_THREADS \
    -n - \
    ${BAM}/${NAME}_*_PROseq_dedup_QC_end.sort.bam | \
  samtools sort \
    -@ $OMP_NUM_THREADS | \
  samtools view \
    -@ $OMP_NUM_THREADS \
    -f 3 -F 2316 \
    -q ${qscore} \
    -o ${mgsam_p}

  # TEXp
  ##samtools merge \
  ##  -@ $OMP_NUM_THREADS \
  ##  -n - \
  ##  ${BAM}/${NAME}_*_TEXp_dedup_QC_end.sort.bam | \
  ##samtools sort \
  ##  -@ $OMP_NUM_THREADS | \
  ##samtools view \
  ##  -@ $OMP_NUM_THREADS \
  ##  -f 3 -F 2316 \
  ##  -q ${qscore} \
  ##  -o ${mgsam_p}

  # TEXm as minus
  samtools merge \
    -@ $OMP_NUM_THREADS \
    -n - \
    ${BAM}/${NAME}_*_TEXm_dedup_QC_end.sort.bam | \
  samtools sort \
    -@ $OMP_NUM_THREADS | \
  samtools view \
    -@ $OMP_NUM_THREADS \
    -f 3 -F 2316 \
    -q ${qscore} \
    -o ${mgsam_m}

else

  echo "${NAME} -- existing merged SAM found"

fi

if [[ ! -f ${CLIST} ]]; then

  echo "${NAME} -- getting mapped contig list"

  awk -F$'\t' '{print $3}' ${mgsam_p} | \
    sort | uniq -c | sed "s/^ *//g" \
    > ${NAME}_readcount_PLUS.txt
  awk -F$'\t' '{print $3}' ${mgsam_m} | \
    sort | uniq -c | sed "s/^ *//g" \
    > ${NAME}_readcount_MINUS.txt

  echo "${NAME} -- filtering contigs by length:${length}, reads:${nreads}"

  join ${NAME}_readcount_PLUS.txt ${NAME}_readcount_MINUS.txt -j 2 | \
    awk -F_ -v len=${length} '$4 >= len' | \
    awk -F' ' -v cut=${nreads} '$2 >= cut && $3 >= cut {print $1}' \
    > ${CLIST}

  rm ${NAME}_readcount_PLUS.txt
  rm ${NAME}_readcount_MINUS.txt

else

  echo "${NAME} -- existing filtered contig list found"

fi

ncon=$(wc -l < ${CLIST})
echo "${NAME} -- ${ncon} contigs passing mapped read threshold"

source ${HOME}/miniconda3/bin/activate
conda activate tssar

while read contig; do

  bcon=$(echo ${contig} | sed "s/_length.*//g")
  mkdir -p ${OUT}/${bcon}/tmp
  cd ${OUT}/${bcon}

  mcon=$(grep -Pn "^${contig}$" ${CLIST} | cut -f1 -d:)
  echo "${NAME} -- running TSSAR on contig ${mcon} / ${ncon}"

  grep -P "\t${contig}\t" ${mgsam_p} > ${contig}_PLUS.sam
  grep -P "\t${contig}\t" ${mgsam_m} > ${contig}_MINUS.sam

  # grep -A1 works with unwrapped fasta only
  grep -A1 ">${contig}" ${SEQ} > ${contig}.fasta

  TSSAR \
    --libP ${contig}_PLUS.sam \
    --libM ${contig}_MINUS.sam \
    --fasta ${contig}.fasta \
    --minPeak 3 \
    --winSize 1000 \
    --range 3 \
    --clean \
    --tmpdir ${OUT}/${bcon}/tmp \
    >> ${OUT}/${NAME}_TSS_${length}nt_${nreads}reads.bed

  rm ${contig}_PLUS.sam
  rm ${contig}_MINUS.sam
  rm ${contig}.fasta
  rm -r ${OUT}/${bcon}/tmp

done < ${CLIST}

echo "${NAME} -- TSSAR done! Removing sam files"

rm ${mgsam_p}
rm ${mgsam_m}

