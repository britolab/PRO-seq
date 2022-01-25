#$ -S /bin/bash
#$ -N groupsums
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq2/log/groupsum_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq2/log/groupsum_$JOB_ID.err
#$ -wd /workdir/users/acv46/ben_May2021
#$ -l h_vmem=80G
#$ -q long.q@cbsubrito
#$ -t 1-2

WRK=/workdir/users/acv46/stool_PROSeq2/deseq/bam/pileup
export OMP_NUM_THREADS=8

DESIGN_FILE=$WRK/files.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

sample=$(echo $NAME | awk -F_ '{print $1}')

cd $WRK

echo "${sample} -- PRO-seq 3p"

grep "3,PRO-seq" $NAME | \
	cut -d, -f4 | \
	sed "s/-//g" | \
	awk '{s+=$1}END{print s}' | \
	sed "s/^/PRO-seq\t3\t/g" \
	>> $WRK/${sample}_groupsums.txt

echo "${sample} -- PRO-seq 5p"

grep "5,PRO-seq" $NAME | \
        cut -d, -f4 | \
        sed "s/-//g" | \
        awk '{s+=$1}END{print s}' | \
        sed "s/^/PRO-seq\t5\t/g" \
        >> $WRK/${sample}_groupsums.txt

echo "${sample} -- PRO-seq full"

grep "full,PRO-seq" $NAME | \
        cut -d, -f4 | \
        sed "s/-//g" | \
        awk '{s+=$1}END{print s}' | \
        sed "s/^/PRO-seq\tfull\t/g" \
        >> $WRK/${sample}_groupsums.txt

echo "${sample} -- RNA-seq full"

grep "full,RNA-seq" $NAME | \
        cut -d, -f4 | \
        sed "s/-//g" | \
        awk '{s+=$1}END{print s}' | \
        sed "s/^/RNA-seq\tfull\t/g" \
        >> $WRK/${sample}_groupsums.txt

echo "${sample} -- done!"
