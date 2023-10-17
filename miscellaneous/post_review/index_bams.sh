#$ -S /bin/bash
#$ -N sort_index_bam
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq3/log/sortindexbams_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq3/log/sortindexbams_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq3
#$ -l h_vmem=20G
#$ -q long.q@cbsubrito
#$ -t 1-22
#$ -tc 4
#$ -pe parenv 4

# sort bam files and create .bai indices
# needed for IGV in VNC Viewer

export OMP_NUM_THREADS=4
export PATH=/programs/samtools-1.15.1-r/bin:$PATH
WRK=/workdir/users/acv46/stool_PROSeq3/transcriptomes

List=${WRK}/bams.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" ${List})
        NAME=`basename "$DESIGN"`

samp=$(echo "${NAME}" | awk -F_ '{print $1}')
short=$(echo "${NAME}" | sed "s/_dedup_QC_end.sort.bam//g")

BAM=${WRK}/${samp}/${NAME}
OUT=${WRK}/${samp}/${short}_resort.bam

cd ${WRK}/${samp}

samtools sort \
  -@ ${OMP_NUM_THREADS} \
  --write-index \
  -o ${OUT} \
  ${BAM}

mkdir -p ${WRK}/${samp}/old_bam
mv ${BAM} ${WRK}/${samp}/old_bam
