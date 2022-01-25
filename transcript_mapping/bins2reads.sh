#$ -S /bin/bash
#$ -N bins2reads
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/bins2reads_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/bins2reads_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=50G
#$ -q long.q@cbsubrito
#$ -t 1-2

# simple script to count the number of reads mapping to specific bins

WRK=/workdir/users/acv46/stool_PROSeq2/deseq/grid
OUT=$WRK/bins2reads

DESIGN_FILE=$WRK/jobs.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

base=$(echo ${NAME} | awk -F_ '{print $1}')

cd $OUT
awk 'NR>1 {print $1}' $WRK/${NAME}_out/${NAME}.interleaved.GRiD.txt \
	> ${NAME}.gridbins.txt

# create text files to index
BAM=/workdir/users/acv46/stool_PROSeq2/deseq/bam
samtools view -f 3 -F 2828 ${BAM}/${base}_RNA1.allcontigs.bam | awk '{print $3}' > $OUT/${NAME}_RNA1.tmp
samtools view -f 3 -F 2828 ${BAM}/${base}_RNA2.allcontigs.bam | awk '{print $3}' > $OUT/${NAME}_RNA2.tmp
samtools view -f 3 -F 2828 ${BAM}/${base}_PRO1.allcontigs.bam | awk '{print $3}' > $OUT/${NAME}_PRO1.tmp
samtools view -f 3 -F 2828 ${BAM}/${base}_PRO2.allcontigs.bam | awk '{print $3}' > $OUT/${NAME}_PRO2.tmp

contigs=/workdir/users/acv46/mgmAssembly/${NAME}_CAB/DASTool/${NAME}*DASTool_scaffolds2bin.txt

while read bin; do

	grep -P "\t${bin}$" ${contigs} | \
	awk -F$'_cov' '{print $1}' \
	> ${NAME}_${bin}_list.tmp

	while read con; do

		rna1=$(grep -w "${con}" $OUT/${NAME}_RNA1.tmp | wc -l)
		rna2=$(grep -w "${con}" $OUT/${NAME}_RNA2.tmp | wc -l)
		pro1=$(grep -w "${con}" $OUT/${NAME}_PRO1.tmp | wc -l)
		pro2=$(grep -w "${con}" $OUT/${NAME}_PRO2.tmp | wc -l)

		echo -e "${NAME}\t${bin}\t${con}\t${rna1}\t${rna2}\t${pro1}\t${pro2}" \
		>> ${NAME}_bin2reads.txt

	done < ${NAME}_${bin}_list.tmp
	rm ${NAME}_${bin}_list.tmp

done < ${NAME}.gridbins.txt

rm $OUT/${NAME}_RNA1.tmp
rm $OUT/${NAME}_RNA2.tmp
rm $OUT/${NAME}_PRO1.tmp
rm $OUT/${NAME}_PRO2.tmp
