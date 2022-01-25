#$ -S /bin/bash
#$ -N tRNA
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/tRNA_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/tRNA_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=150G
#$ -q long.q@cbsubrito
#$ -t 1

SCRIPT=/workdir/users/acv46/stool_PROSeq2/scripts
PATH=/programs/R-4.1.2/bin:$PATH

Rscript $SCRIPT/tRNA_coverage.R
