#$ -S /bin/bash
#$ -N antismash
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/antismash_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/antismash_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=100G
#$ -q long.q@cbsubrito
#$ -t 1-2

# source /home/acv46/miniconda3/bin/activate
# conda create -n antismash
# conda activate antismash
# conda install -c bioconda antismash
# download-antismash-databases
# conda install natsort

WRK=/workdir/users/acv46/stool_PROSeq2
ANTI=/programs/bin/run_antismash
export OMP_NUM_THREADS=8

DESIGN_FILE=$WRK/antismash/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

samp=$(echo $NAME | awk -F_ '{print $1}')
gbk=$WRK/deseq/prokka/${samp}_allcontigs/${samp}_allcontigs.gbk

OUT=$WRK/antismash/${NAME}_out
mkdir $OUT
cd $OUT
GFF=$WRK/deseq/prokka/${samp}_allcontigs/${samp}_allcontigs.gff

# note: output dir cannot contain other files

source /home/acv46/miniconda3/bin/activate
conda activate antismash

antismash \
	--taxon bacteria \
	-c $OMP_NUM_THREADS \
	--output-dir $OUT \
	--output-basename $samp \
	--genefinding-tool prodigal-m \
	$gbk

conda deactivate
