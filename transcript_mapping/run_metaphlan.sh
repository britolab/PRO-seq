#$ -S /bin/bash
#$ -N metaphlan
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/proseq_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/proseq_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=80G
#$ -q long.q@cbsubrito
#$ -t 1-8

WRK=/workdir/users/acv46/stool_PROSeq2
JOB=$WRK/scripts/samples.txt

LIST=$JOB
  DESIGN=$(sed -n "${SGE_TASK_ID}p" $LIST)
  NAME=`basename "$DESIGN"`

echo "${NAME} -- starting metaphlan"

export PATH=/programs/MetaPhlAn-2.0:/programs/MetaPhlAn-2.0/utils:$PATH
export OMP_NUM_THREADS=8

FASTQ1=$WRK/fastq/clean/${NAME}_R1.fastq.gz
FASTQ2=$WRK/fastq/clean/${NAME}_R2.fastq.gz

MPHL=$WRK/out/${NAME}/metaphlan
mkdir -p $MPHL
cd $MPHL

metaphlan2.py \
        --input_type fastq \
        --bowtie2db /workdir/users/fnn3/references/db_v20/mpa_v20_m200 \
        --bowtie2out ${NAME}_mgm.bowtie2.bz2 \
        --nproc $OMP_NUM_THREADS \
        $FASTQ1,$FASTQ2 \
        ${NAME}_profile.txt

echo "${NAME} -- metaphlan finished"

