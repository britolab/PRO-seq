#$ -S /bin/bash
#$ -N hicpro
#$ -V
#$ -t 1
#$ -e /workdir/users/acv46/hidin-seq/log/hicpro_$JOB_ID.err
#$ -o /workdir/users/acv46/hidin-seq/log/hicpro_$JOB_ID.out
#$ -wd /workdir/users/acv46/hidin-seq
#$ -l h_vmem=30G
#$ -q long.q@cbsubrito

#This script runs hicpro setup.
WRK=/workdir/users/acv46/hidin-seq
#DESIGN_FILE=$WRK/HicDesign.txt
#        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
#        NAME=`basename "$DESIGN"`

SCF=$WRK/US3_10/bins/bin7_plas29482_ref.fasta
SCFNAME=US3_10_bin7_plas29482
cd $WRK/US3_10

#export PYTHONPATH=/home/britolab/acv46/miniconda2/lib/python2.7/site-packages:$PYTHONPATH
#export PATH=/programs/HiC-Pro_2.7.9/bin:$PATH

export PYTHONPATH=/programs/HiC-Pro_2.7.9/lib64/python2.7/site-packages:$PYTHONPATH
export PATH=/programs/HiC-Pro_2.7.9/bin:$PATH

#build bowtie2 index
bowtie2-build $SCF $SCFNAME

#get chromosome name & size file
#python ~/agk/CDC/scripts/reference_genome_sizes.py $SCF ${NAME}_scaffold.sizes

#get bed file of restriction enzyme sites
python2.7 /programs/HiC-Pro_2.7.9/bin/utils/digest_genome.py \
-r ^GATC \
-o bin7_plas29482_DpnII.bed \
$SCF

#OUTPUT=$WRK/hicpro/${NAME}_rawdata/${NAME}
#if [ ! -d $OUTPUT ]; then mkdir -p $OUTPUT; fi
#R1=/workdir/data/CDC/hic/merged/${NAME}hic.1.fastq
#R2=/workdir/data/CDC/hic/merged/${NAME}hic.2.fastq
#cp $R1 ${OUTPUT}/${NAME}_R1.fastq
#cp $R2 ${OUTPUT}/${NAME}_R2.fastq
#mkdir $WRK/hicpro/${NAME}_rawdata
#/programs/HiC-Pro_2.7.9/bin/utils/split_reads.py --results_folder $WRK/hicpro/${NAME}_rawdata/${NAME} --nreads 500000 $OUTPUT/${NAME}_R1.fastq
#/programs/HiC-Pro_2.7.9/bin/utils/split_reads.py --results_folder $WRK/hicpro/${NAME}_rawdata/${NAME} --nreads 500000 $OUTPUT/${NAME}_R2.fastq

#cp /workdir/users/agk85/CDC2/hicpro/configs/config-hicpro.txt /workdir/users/agk85/CDC2/hicpro/configs/${NAME}_config-hicpro.txt
#sed -i "s/B331-2/$NAME/g" "/workdir/users/agk85/CDC2/hicpro/configs/${NAME}_config-hicpro.txt"

