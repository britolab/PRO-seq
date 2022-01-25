#$ -S /bin/bash
#$ -N vibrant2propagate
#$ -V
#$ -o /workdir/users/acv46/ben_May2021/log/vibrant2propagate_$JOB_ID.out
#$ -e /workdir/users/acv46/ben_May2021/log/vibrant2propagate_$JOB_ID.err
#$ -wd /workdir/users/acv46/ben_May2021
#$ -l h_vmem=80G
#$ -q long.q@cbsubrito2
#$ -t 1-6
#$ -m ea
#$ -M acv46@cornell.edu

# this script first runs vibrant on metagenomic contigs
# then, runs propagate on prophage coordinates

# both programs installed via conda
# see /home/acv46/propagate/_README
# and /home/acv46/vibrant/_README

WRK=/workdir/users/acv46/ben_May2021
export OMP_NUM_THREADS=8
export PATH=/programs/samtools-1.11/bin:$PATH
export PATH=/programs/bowtie2-2.4.3-linux-x86_64:$PATH

vout=$WRK/phage/vibrant
VIB=/home/acv46/vibrant/VIBRANT/VIBRANT_run.py
PRO=/home/acv46/propagate/PropagAtE-main/PropagAtE_run.py
fastq=/home/britolab/data/2021_04_03_nextseq_plaque_metagenome

source /home/acv46/miniconda3/bin/activate
conda activate vibrant

DESIGN_FILE=$WRK/vamb/small_NN/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

mkdir -p $vout
cd $vout

echo "${NAME} -- starting vibrant"

python3 $VIB \
	-i $WRK/vamb/small_NN/bins_merged_samples/${NAME}_merged.fna \
	-f nucl \
	-folder $vout/${NAME} \
	-t $OMP_NUM_THREADS \
	-l 5000 \
	-o 10

echo "${NAME} -- finished vibrant"

conda activate propagate

echo "${NAME} -- starting propagate"

cd $vout/${NAME}

index=$WRK/vamb/vamb_rename.txt
samp=$(grep "${NAME}" $index | awk '{print $1}')

python3 $PRO \
	-v VIBRANT_${NAME}_merged/VIBRANT_results_${NAME}_merged/VIBRANT_integrated_prophage_coordinates_${NAME}_merged.tsv \
	-r $fastq/${samp}.clean_1.fastq.gz $fastq/${samp}.clean_2.fastq.gz \
	-f $WRK/vamb/small_NN/bins_merged_samples/${NAME}_merged.fna \
	-t $OMP_NUM_THREADS \
	-o ${NAME}_propagate_out.tsv \
	-clean

conda deactivate

echo "${NAME} -- finished propagate, script done"
