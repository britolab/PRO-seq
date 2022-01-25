#$ -S /bin/bash
#$ -N phasefinder
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/phasefinder_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/phasefinder_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=100G
#$ -q long.q@cbsubrito
#$ -t 1-2

# runs the PhaseFinder script, using RNA-seq and PRO-seq inputs

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

WRK=/workdir/users/acv46/stool_PROSeq2/invertons/out_mismatch0
pfind=/home/britolab/acv46/PhaseFinder/PhaseFinder.py

DESIGN_FILE=$WRK/../jobs.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

MGM=/workdir/users/acv46/mgmAssembly/${NAME}_CAB/metaSPAdes/contigs.fasta
FQ1R1=/workdir/users/acv46/stool_PROSeq2/fastq/clean/${NAME}_RNA_rep1_R1.fastq.gz
FQ2R1=/workdir/users/acv46/stool_PROSeq2/fastq/clean/${NAME}_RNA_rep2_R1.fastq.gz
FQ1R2=/workdir/users/acv46/stool_PROSeq2/fastq/clean/${NAME}_RNA_rep1_R2.fastq.gz
FQ2R2=/workdir/users/acv46/stool_PROSeq2/fastq/clean/${NAME}_RNA_rep2_R2.fastq.gz

FQ1P1=/workdir/users/acv46/stool_PROSeq2/fastq/clean/${NAME}_PRO_rep1_R1.fastq.gz
FQ2P1=/workdir/users/acv46/stool_PROSeq2/fastq/clean/${NAME}_PRO_rep2_R1.fastq.gz
FQ1P2=/workdir/users/acv46/stool_PROSeq2/fastq/clean/${NAME}_PRO_rep1_R2.fastq.gz
FQ2P2=/workdir/users/acv46/stool_PROSeq2/fastq/clean/${NAME}_PRO_rep2_R2.fastq.gz

cd $WRK

# merge replicates and unzip
cat $FQ1R1 $FQ2R1 | gunzip -c > ${NAME}_RNA_R1_merged.fastq
cat $FQ1R2 $FQ2R2 | gunzip -c > ${NAME}_RNA_R2_merged.fastq
cat $FQ1P1 $FQ2P1 | gunzip -c > ${NAME}_PRO_R1_merged.fastq
cat $FQ1P2 $FQ2P2 | gunzip -c > ${NAME}_PRO_R2_merged.fastq

source /home/britolab/acv46/miniconda3/bin/activate
conda activate PhaseFinder
export PATH=/programs/EMBOSS-6.6.0/bin:$PATH
export PATH=/programs/bowtie2-2.3.0:$PATH
export PATH=/programs/bedops-2.4.35/bin:$PATH
export PATH=/programs/bedtools2-2.29.2/bin:$PATH
export OMP_NUM_THREADS=16

# use phase-invertible calls from previous run

echo "${NAME} -- running PhaseFinder ratio with RNA-seq reads"

python $pfind ratio \
	-i $WRK/${NAME}.ID.fasta \
	-1 ${NAME}_RNA_R1_merged.fastq \
	-2 ${NAME}_RNA_R2_merged.fastq \
	-p $OMP_NUM_THREADS \
	-o $WRK/${NAME}_RNA_out

echo "${NAME} -- running PhaseFinder ratio with PRO-seq reads"

python $pfind ratio \
        -i $WRK/${NAME}.ID.fasta \
        -1 ${NAME}_PRO_R1_merged.fastq \
        -2 ${NAME}_PRO_R2_merged.fastq \
        -p $OMP_NUM_THREADS \
        -o $WRK/${NAME}_PRO_out

echo "${NAME} -- PhaseFinder done"

rm ${NAME}_RNA_R1_merged.fastq
rm ${NAME}_RNA_R2_merged.fastq
rm ${NAME}_PRO_R1_merged.fastq
rm ${NAME}_PRO_R2_merged.fastq
