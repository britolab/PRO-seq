#$ -S /bin/bash
#$ -N prodigal
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/prodigal_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/prodigal_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=50G
#$ -q long.q@cbsubrito
#$ -t 1-4

WRK=/workdir/users/acv46/stool_PROSeq2/deseq
OUT=$WRK/prodigal
export PATH=/programs/prodigal-2.6.3:$PATH

DESIGN_FILE=$WRK/jobs.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

cd $OUT

# simplify fasta header
fasta=$(grep "${NAME}" $WRK/ref.txt | awk '{print $2}')
sed 's/_cov.*//g' $fasta > ${NAME}.short.fasta

prodigal \
	-i ${NAME}.short.fasta \
	-p meta \
	-o ${NAME}.gff \
	-f gff \
	-g 11 \
	-d ${NAME}_mRNA.fasta

