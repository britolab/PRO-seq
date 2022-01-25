#$ -S /bin/bash
#$ -N parse_antismash
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/parse_antismash_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/parse_antismash_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=5G
#$ -q long.q@cbsubrito
#$ -t 1-2

WRK=/workdir/users/acv46/stool_PROSeq2/antismash
export OMP_NUM_THREADS=8

DESIGN_FILE=$WRK/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

base=$(echo ${NAME} | awk -F_ '{print $1}')
list=$(ls $WRK/${NAME}_out | grep "region001.gbk")

cd $WRK/${NAME}_out

while read gbk; do

	start=$(grep "Orig. start" $gbk | tr -d -c 0-9)
	end=$(grep "Orig. end" $gbk | tr -d -c 0-9)
	contig=$(grep "Original ID" $gbk | awk -F:: '{print $2}' | sed "s/ //g")
	genelist=$(grep -P "^ *gene " $gbk)

	while read gene; do

		gstart=$(echo $gene | sed "s/gene //g" | sed "s/complement(//g" | awk -F. '{print $1}')
		gend=$(echo $gene | sed "s/gene //g" | sed "s/complement(//g" | awk -F. '{print $3}' | sed "s/)//g")
		gname=$(grep -A2 "$gene" $gbk | grep "/gene" | sed 's|/gene="||g' | sed 's/"//g' | tr -d [[:blank:]] | awk -F_ '{print $1}')

		if [ -z "$gname" ]; then
			gname=NA
		fi

		echo $gene | grep -q "complement"
		if [ $? -eq 0 ]; then
			strand=-1
		else
			strand=1
		fi

		echo -e "$gbk\t$start\t$end\t$contig\t$gstart\t$gend\t$gname\t$strand" \
		>> ${base}_temp.txt

	done <<< "$genelist"

done <<< "$list"

grep -v "/gene" ${base}_temp.txt > ${base}_antismash_summary.txt
rm ${base}_temp.txt
