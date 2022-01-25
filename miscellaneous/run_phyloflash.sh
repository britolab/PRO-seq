#$ -S /bin/bash
#$ -N phyloflash
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/phyloflash_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/phyloflash_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=80G
#$ -q long.q@cbsubrito
#$ -t 1-2

# phyloflash is a pipeline for reconstruction of small subunit ribosomal RNAs

# installation
# https://github.com/HRGV/phyloFlash
## source /home/britolab/acv46/miniconda3/bin/activate
## conda create -n phyloflash
## conda activate phyloflash
## cd $HOME
## wget https://github.com/HRGV/phyloFlash/archive/pf3.4.tar.gz
## tar -xzf pf3.4.tar.gz
## cd phyloFlash-pf3.4
## ./phyloFlash.pl -check_env

WRK=/workdir/users/acv46/stool_PROSeq2
PFL=$HOME/phyloFlash-pf3.4
OUT=$WRK/phyloflash
DB=$PFL/138.1
export OMP_NUM_THREADS=16
export PATH=/programs/spades/bin:$PATH
export PATH=/programs/bbmap-38.90:$PATH
export PATH=/programs/bbmap-38.90/reformat.sh:$PATH
export PATH=/programs/mafft/bin:$PATH
export PATH=/programs/vsearch-2.15.0/bin:$PATH


DESIGN_FILE=$OUT/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

base=$(echo $NAME | cut -f1 -d "_")

mkdir -p $OUT/${NAME}
cd $OUT/${NAME}

source /home/britolab/acv46/miniconda3/bin/activate
conda activate phyloflash

# run spades and EMIRGE for 16S rRNA sequence reconstruction

FQ1=/workdir/users/acv46/stool_PROSeq2/fastq/clean/${NAME}_mgm_R1.fastq
FQ2=/workdir/users/acv46/stool_PROSeq2/fastq/clean/${NAME}_mgm_R2.fastq

# need to format db correctly to use sortmerna

$PFL/phyloFlash.pl \
	-lib $base \
	-everything \
	-read1 $FQ1 \
	-read2 $FQ2 \
	-CPUs $OMP_NUM_THREADS \
	-dbhome $DB
	##-sortmerna

