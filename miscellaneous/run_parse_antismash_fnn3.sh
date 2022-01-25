#$ -S /bin/bash
#$ -N parse_antismash
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/parse_antismash_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/parse_antismash_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=5G
#$ -q long.q@cbsubrito
#$ -t 1-2

#This script will parse the antismash output genbank files
#Python script is from Nielsen et al 2017

WRK=/workdir/users/acv46/stool_PROSeq2/antismash
parse=/workdir/users/acv46/stool_PROSeq2/scripts/antismash_gbkParse.py

DESIGN_FILE=$WRK/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

base=$(echo ${NAME} | awk -F_ '{print $1}')
cd $WRK/${NAME}_out
python2 $parse ${base}.gbk
