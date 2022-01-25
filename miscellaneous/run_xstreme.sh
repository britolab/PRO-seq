#$ -S /bin/bash
#$ -N xstreme
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/xstreme_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/xstreme_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=50G
#$ -q long.q@cbsubrito
#$ -t 1-2

# MEME suite version 5.4.1
# see installation guide at
# https://meme-suite.org/meme/doc/install.html?man_type=web

WRK=/workdir/users/acv46/stool_PROSeq2
OUT=$WRK/dbcan2
MEME=/home/acv46/meme/bin
export OMP_NUM_THREADS=4

DESIGN_FILE=$OUT/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

$MEME/xstreme --p \
	--oc \
	--dna \
	--minw 6 \
	--maxw 20 \
	--verbosity 3
