#$ -S /bin/bash
#$ -N NNPP2.2
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/NNPP2.2_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/NNPP2.2_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=20G
#$ -q long.q@cbsubrito
#$ -t 1-2

export OMP_NUM_THREADS=16
WRK=/workdir/users/acv46/stool_PROSeq2
IN=$WRK/deseq/ref
OUT=$WRK/promoters
NNPP22=/home/acv46/MetaRon-master/NNPP2.2/bin/fa2TDNNpred-PRO.linux

List=$OUT/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $List)
        NAME=`basename "$DESIGN"`

cd $OUT

awk -v name="${NAME}.fa" '/^>/{s=++d} {print > s "_" name}' \
	${IN}/${NAME}_allcontigs.short_1000.fasta

for file in *_${NAME}.fa; do

	head=$(grep ">" ${file} | sed "s/>//g")
	${NNPP22} -r -t 0.9 $file | \
	sed "s/^Hit/${head} - Hit/g" | \
	grep -A2 "position:" | \
	tr '\n' ',' | sed "s/---*/\n/g" | \
	sed "s/,NODE/NODE/g" | grep "NODE" | \
	awk -F$" " '{print $1,$4, $6, $8, $12, $14, $15}' | \
	sed "s/://g" | sed "s/,//g" | sed "s/Prediction//g" \
	>> ${NAME}_promoters.txt

	rm ${file}

done
