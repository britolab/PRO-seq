#$ -S /bin/bash
#$ -N minced
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq3/log/minced_Ec_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq3/log/minced_Ec_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq3
#$ -l h_vmem=50G
#$ -q long.q@cbsubrito
#$ -t 1
#$ -tc 1
#$ -pe parenv 4

# runs minced
# installed via makefile

WRK=/workdir/users/acv46/stool_PROSeq3
#BLASTDB=$WRK/phage/virsorter
MIN=/home/acv46/minced/minced
OUT=$WRK/figures/Ecoli

export OMP_NUM_THREADS=4

mkdir -p $OUT/minced
cd $OUT/minced

$MIN -gffFull -spacers \
	-minNR 3\
	-minRL 20 \
	-maxRL 50 \
	-minSL 22 \
	-maxSL 55 \
	$OUT/GCF_000005845.2_ASM584v2_genomic.fna \
	ASM584v2_minced.txt \
	ASM584v2_minced.gff

