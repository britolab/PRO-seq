#$ -S /bin/bash
#$ -N prokka
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/prokka_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/prokka_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=125G
#$ -q long.q@cbsubrito
#$ -t 1-2

## in miniconda3

# conda install -c conda-forge -c bioconda -c defaults prokka

WRK=/workdir/users/acv46/stool_PROSeq2/out3
export OMP_NUM_THREADS=12

DESIGN_FILE=$WRK/jobs.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

OUT=$WRK/${NAME}_ref
mkdir -p $OUT

# run prokka for highly fragmented genomes (--metagenome) and include ncRNAs (--rfam)

source $HOME/miniconda3/bin/activate

cd $OUT

# simplify fasta header
# fasta=$(grep "${NAME}" $WRK/ref.txt | awk '{print $2}')
# sed 's/_cov.*//g' $fasta > ${NAME}.short.fasta

prokka \
	--cpus $OMP_NUM_THREADS \
	--metagenome \
	--rfam \
	--gffver 3 \
	--gcode 11 \
	--addgenes \
	--prefix ${NAME}_bins_plasmids \
	${NAME}_bins_plasmids.fa

