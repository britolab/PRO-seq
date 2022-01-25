#$ -S /bin/bash
#$ -N phasefinder
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/phasefinder_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/phasefinder_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=100G
#$ -q long.q@cbsubrito
#$ -t 1-2

# create PhaseFinder environment:

## source /home/britolab/acv46/miniconda3/bin/activate
## conda create --name PhaseFinder python=2.7
## conda activate PhaseFinder
## conda install biopython
## conda install pandas

# Run PhaseFinder:

## source /home/britolab/acv46/miniconda3/bin/activate
## conda activate PhaseFinder
## export PATH=/programs/EMBOSS-6.6.0/bin:$PATH
## export PATH=/programs/bowtie2-2.3.0:$PATH
## export PATH=/programs/bedops-2.4.35/bin:$PATH
## export PATH=/programs/bedtools2-2.29.2/bin:$PATH
## python PhaseFinder.py -h

WRK=/workdir/users/acv46/stool_PROSeq2/invertons
pfind=/home/britolab/acv46/PhaseFinder/PhaseFinder.py

DESIGN_FILE=$WRK/jobs.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

MGM=/workdir/users/acv46/mgmAssembly/${NAME}_CAB/metaSPAdes/contigs.fasta
FQ1=/workdir/users/acv46/stool_PROSeq2/fastq/clean/${NAME}_mgm_R1.fastq
FQ2=/workdir/users/acv46/stool_PROSeq2/fastq/clean/${NAME}_mgm_R2.fastq

source /home/britolab/acv46/miniconda3/bin/activate
conda activate PhaseFinder
export PATH=/programs/EMBOSS-6.6.0/bin:$PATH
export PATH=/programs/bowtie2-2.3.0:$PATH
export PATH=/programs/bedops-2.4.35/bin:$PATH
export PATH=/programs/bedtools2-2.29.2/bin:$PATH
export OMP_NUM_THREADS=16

echo "${NAME} -- running PhaseFinder locate"

python $pfind locate \
	-f $MGM \
	-t $WRK/${NAME}.einverted.tab \
	-g 15 85 \
	-p \
	-m 1

echo "${NAME} -- running PhaseFinder create"

python $pfind create \
	-f $MGM \
	-t $WRK/${NAME}.einverted.tab \
	-s 1000 \
	-i $WRK/${NAME}.ID.fasta

echo "${NAME} -- running PhaseFinder ratio"

python $pfind ratio \
	-i $WRK/${NAME}.ID.fasta \
	-1 $FQ1 \
	-2 $FQ2 \
	-p $OMP_NUM_THREADS \
	-o $WRK/${NAME}_out

echo "${NAME} -- PhaseFinder done"
