#$ -S /bin/bash
#$ -N GRiD
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/GRiD_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/GRiD_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=100G
#$ -q long.q@cbsubrito
#$ -t 1-2

# conda installation
## source /home/britolab/acv46/miniconda3/bin/activate
## conda config --add channels defaults
## conda config --add channels bioconda
## conda config --add channels conda-forge
## conda create --name GRiD
## conda activate GRiD
## conda install grid=1.3

# custom DASTool bin databases created with the following command
# cd /workdir/users/acv46/stool_PROSeq2/deseq/grid/US3_3Nov2020
# update_database -d . -g /workdir/users/acv46/mgmAssembly/US3_3Nov2020_CAB/DASTool/US3_3Nov2020_DASTool_bins -p US3_3Nov2020

WRK=/workdir/users/acv46/stool_PROSeq2/deseq/grid

DESIGN_FILE=$WRK/jobs.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

# concatenate paired-end reads (GRiP only takes single-end reads)

echo "${NAME} -- concatenating PE fastqs"

fq=/workdir/users/acv46/mgmAssembly/${NAME}_CAB/fastq
mkdir -p ${fq}/interleaved
/programs/bbmap-38.86/reformat.sh \
	in1=${fq}/${NAME}.clean_1.fastq \
	in2=${fq}/${NAME}.clean_2.fastq \
	out=${fq}/interleaved/${NAME}.interleaved.fastq

# run grip multiplex

echo "${NAME} -- starting GRiD using database at $WRK/db/${NAME}"

source /home/britolab/acv46/miniconda3/bin/activate
conda activate GRiD

export OMP_NUM_THREADS=16
export PATH=/programs/bowtie2-2.3.0:$PATH
export PATH=/programs/seqtk:$PATH
export PATH=/programs/R-3.4.1/bin:$PATH
export PATH=/programs/samtools-1.6/bin:$PATH
export PATH=/programs/bedtools2-2.26.0/bin:$PATH
#export PATH=/programs/bamtools-2.2.3/bin:$PATH
export PATH=/programs/ncbi-blast-2.9.0+/bin:$PATH

grid multiplex \
	-r ${fq}/interleaved \
	-d $WRK/db/${NAME} \
	-p \
	-c 0.2 \
	-o ${WRK}/${NAME}_out \
	-n ${OMP_NUM_THREADS}

# removed interleaved fastq to save space

echo "${NAME} -- GRiD finished, removing interleaved fastq"

rm -r ${fq}/interleaved

echo "${NAME} -- GRiD pipeline DONE"
