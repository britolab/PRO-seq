#$ -S /bin/bash
#$ -N kraken2
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq3/log/kraken2_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq3/log/kraken2_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq3
#$ -l h_vmem=65G
#$ -q long.q@cbsubrito
#$ -t 1-10
#$ -tc 3
#$ -M acv46@cornell.edu
#$ -pe parenv 8

# kraken2 installed via conda
## source /home/acv46/miniconda3/bin/activate
## conda create --name kraken2
## conda activate kraken2
## conda install -c bioconda kraken2

WRK=/workdir/users/acv46/stool_PROSeq3/kraken
export OMP_NUM_THREADS=8
# see /workdir/refdbs/kraken/HumGut_DB/_README for database install details
kdb=/workdir/refdbs/kraken2db
kout=$WRK/kraken2_old_bam

source /home/acv46/miniconda3/bin/activate
conda activate kraken2

DESIGN_FILE=$WRK/old_names.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

OUT=$kout/${NAME}
mkdir -p $OUT

bam=$(grep "${NAME}" $WRK/old_names2bams.txt | cut -f2)

echo "${NAME} -- starting kraken2"

cd $OUT

samtools fastq ${bam} > ${NAME}_converted.fastq

kraken2 \
	--report ${NAME}.report.txt \
	--db ${kdb} \
	--threads ${OMP_NUM_THREADS} \
	--output ${NAME}.out.txt \
	${NAME}_converted.fastq

echo "${NAME} -- finished kraken2"

conda deactivate
