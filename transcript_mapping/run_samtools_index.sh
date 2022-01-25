#$ -S /bin/bash
#$ -N get_CRISPR_insert
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/samtools_index_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/samtools_index_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=50G
#$ -q long.q@cbsubrito
#$ -t 1-4

WRK=/workdir/users/acv46/stool_PROSeq2/deseq/bam/pileup
export OMP_NUM_THREADS=12
export PATH=/programs/samtools-1.11/bin:$PATH

cd $WRK

DESIGN_FILE=$WRK/bams.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

echo "${NAME} -- starting index"

samtools index -@ $OMP_NUM_THREADS ${NAME}

echo "${NAME} -- finished index"
