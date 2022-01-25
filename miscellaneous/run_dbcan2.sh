#$ -S /bin/bash
#$ -N dbcan2
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/dbcan2_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/dbcan2_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=50G
#$ -q long.q@cbsubrito
#$ -t 1-2

# conda installation
# https://github.com/linnabrown/run_dbcan
## source /home/britolab/acv46/miniconda3/bin/activate
## conda create -n dbcan2.0.11 python=3.8 diamond hmmer prodigal -c conda-forge -c bioconda
## conda activate dbcan2.0.11
## pip install run-dbcan==2.0.11

WRK=/workdir/users/acv46/stool_PROSeq2
OUT=$WRK/dbcan2
export OMP_NUM_THREADS=4

DESIGN_FILE=$OUT/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

base=$(echo $NAME | cut -f1 -d "_")

#fna=/home/acv46/dbcan2/EscheriaColiK12MG1655.fna

faa=$WRK/deseq/prokka/${base}_allcontigs/${base}_allcontigs.faa
gff=$WRK/deseq/prokka/${base}_allcontigs/${base}_allcontigs.gff

mkdir $OUT/${NAME}
cd $OUT/${NAME}

source /home/britolab/acv46/miniconda3/bin/activate
conda activate dbcan2.0.11

# protein setting using existing prokka annotations

run_dbcan.py $faa protein \
        -c $gff \
	--out_dir $OUT/${NAME} \
        --dia_cpu $OMP_NUM_THREADS \
        --hmm_cpu $OMP_NUM_THREADS \
        --hotpep_cpu $OMP_NUM_THREADS \
        --tf_cpu $OMP_NUM_THREADS \
        --db_dir /home/acv46/dbcan2/db
