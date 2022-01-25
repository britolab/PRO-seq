#$ -S /bin/bash
#$ -N pasific
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/pasific_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/pasific_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2/pasific
#$ -l h_vmem=30G
#$ -q long.q@cbsubrito
#$ -t 1

# dependencies
## perl
## The ViennaRNA Package (https://www.tbi.univie.ac.at/RNA/)
## Rscript
## randomForest R-package (https://cran.r-project.org/web/packages/randomForest/)

export PATH=/programs/ViennaRNA-2.2.4/bin:$PATH
export PATH=/programs/R-4.0.0/bin:$PATH
PAS=/home/acv46/pasific/PASIFIC

out=/workdir/users/acv46/stool_PROSeq2/pasific/
test=US2_allcontigs

perl $PAS/structure_prediction.pl $out $test 1 1 1 1
