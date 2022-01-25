#$ -S /bin/bash
#$ -N vibrant2propagate
#$ -V
#$ -o /workdir/users/acv46/ben_May2021/log/vibrant2propagate_$JOB_ID.out
#$ -e /workdir/users/acv46/ben_May2021/log/vibrant2propagate_$JOB_ID.err
#$ -wd /workdir/users/acv46/ben_May2021
#$ -l h_vmem=120G
#$ -q long.q@cbsubrito2
#$ -t 4
#$ -m ea
#$ -M acv46@cornell.edu

# this script first runs vibrant on metagenomic contigs
# then, runs propagate on prophage coordinates

# both programs installed via conda
# see /home/acv46/propagate/_README
# and /home/acv46/vibrant/_README

WRK=/workdir/users/acv46/ben_May2021
export OMP_NUM_THREADS=12
export PATH=/programs/samtools-1.11/bin:$PATH
export PATH=/programs/bowtie2-2.4.3-linux-x86_64:$PATH

vout=$WRK/phage/vibrant
VIB=/home/acv46/vibrant/VIBRANT/VIBRANT_run.py
PRO=/home/acv46/propagate/PropagAtE-main/PropagAtE_run.py
fastq=/home/britolab/data/2021_04_03_nextseq_plaque_metagenome

source /home/acv46/miniconda3/bin/activate
conda activate vibrant

DESIGN_FILE=$WRK/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

mkdir -p $vout
cd $vout

echo "${NAME} -- starting vibrant"

python3 $VIB \
	-i $WRK/${NAME}_CAB/metaSPAdes/contigs.fasta \
	-f nucl \
	-folder $vout/${NAME} \
	-t $OMP_NUM_THREADS \
	-l 5000 \
	-o 10

echo "${NAME} -- finished vibrant"

conda activate propagate

echo "${NAME} -- starting propagate"

cd $vout/${NAME}

python3 $PRO \
	-v VIBRANT_contigs/VIBRANT_results_contigs/VIBRANT_integrated_prophage_coordinates_contigs.tsv \
	-r $fastq/${NAME}.clean_1.fastq.gz $fastq/${NAME}.clean_2.fastq.gz \
	-f $WRK/${NAME}_CAB/metaSPAdes/contigs.fasta \
	-t $OMP_NUM_THREADS \
	-o ${NAME}_propagate_out.tsv \
	-clean

conda deactivate

echo "${NAME} -- finished propagate, script done"
