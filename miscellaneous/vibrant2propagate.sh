#$ -S /bin/bash
#$ -N vibrant2propagate
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq2/log/vibrant2propagate_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq2/log/vibrant2propagate_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=80G
#$ -q long.q@cbsubrito
#$ -t 1-2
#$ -m ea
#$ -M acv46@cornell.edu

# this script first runs vibrant on metagenomic contigs
# then, runs propagate on prophage coordinates

# both programs installed via conda
# see /home/acv46/propagate/_README
# and /home/acv46/vibrant/_README

WRK=/workdir/users/acv46/stool_PROSeq2
export OMP_NUM_THREADS=12
export PATH=/programs/samtools-1.11/bin:$PATH
export PATH=/programs/bowtie2-2.4.3-linux-x86_64:$PATH

vout=$WRK/phage/vibrant
VIB=/home/acv46/vibrant/VIBRANT/VIBRANT_run.py
PRO=/home/acv46/propagate/PropagAtE-main/PropagAtE_run.py
fastq=/workdir/users/acv46/stool_PROSeq2/fastq/clean

source /home/acv46/miniconda3/bin/activate
conda activate vibrant

DESIGN_FILE=$WRK/phage/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

mkdir -p $vout
cd $vout

echo "${NAME} -- starting vibrant"

base=$(echo ${NAME} | awk -F_ '{print $1}')

python3 $VIB \
	-i $WRK/deseq/ref/${base}_allcontigs.short.fasta \
	-f nucl \
	-folder $vout/${base} \
	-t $OMP_NUM_THREADS \
	-l 5000 \
	-o 10

echo "${NAME} -- finished vibrant"

conda activate propagate

echo "${NAME} -- starting propagate"

cd $vout/${base}

python3 $PRO \
	-v VIBRANT_${base}_allcontigs.short/VIBRANT_results_${base}_allcontigs.short/VIBRANT_integrated_prophage_coordinates_*.tsv \
	-r $fastq/${NAME}_mgm_R1.fastq $fastq/${NAME}_mgm_R2.fastq \
	-f $WRK/deseq/ref/${base}_allcontigs.short.fasta \
	-t $OMP_NUM_THREADS \
	-o ${base}_propagate_out.tsv \
	-clean

conda deactivate

echo "${NAME} -- finished propagate, script done"
