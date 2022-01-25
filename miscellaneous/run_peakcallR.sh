#$ -S /bin/bash
#$ -N peakcallR
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/peakcallR_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/peakcallR_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2/peakcallR
#$ -l h_vmem=80G
#$ -q long.q@cbsubrito
#$ -t 1-2

# script to run peakcallR.R

export PATH=/programs/R-4.0.0/bin:$PATH

# the following packages should be installed
## optparse
## tidyverse
## stats
## Biostrings
## pbapply

WRK=/workdir/users/acv46/stool_PROSeq2
OUT=$WRK/peakcallR
PKR=$WRK/scripts/peakcallR.R

DESIGN_FILE=$OUT/jobs.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

Rscript $PKR \
	-o ${NAME}_out \
	-d $OUT \
	-c $OUT/${NAME}_cov_pro_90_5_reduced.csv \
	-b $OUT/${NAME}_bindat.csv \
	-f $WRK/deseq/prokka/${NAME}_allcontigs/${NAME}_allcontigs.fna \
	-e 3 \
	-r 25 \
	-y 10 \
	-z 3 \
	-s 20 \
	-m z
