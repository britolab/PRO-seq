#$ -S /bin/bash
#$ -N prevo_bam
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq2/log/prevobam_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq2/log/prevobam_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=50G
#$ -q long.q@cbsubrito
#$ -t 1-5

WRK=/workdir/users/acv46/stool_PROSeq2/anvio/US2_prevo
BINDEX=/workdir/users/acv46/stool_PROSeq2/ref/US2_bins/merged/prevo_bins
FQ=/workdir/users/acv46/stool_PROSeq2/fastq/clean
export OMP_NUM_THREADS=4

DESIGN_FILE=$WRK/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

cd $WRK

new=$(echo $NAME | sed 's/_rep//g')

bwa mem -a -v 2 -t $OMP_NUM_THREADS $BINDEX \
	${FQ}/US2_5Nov2020_${NAME}_R1.fastq.gz \
	${FQ}/US2_5Nov2020_${NAME}_R2.fastq.gz \
	> bam/${new}.unsort.bam

samtools view -u bam/${new}.unsort.bam | \
	samtools sort -@ $OMP_NUM_THREADS -o bam/${new}.bam

rm bam/${new}.unsort.bam

samtools index -@ 4 -b bam/${new}.bam
