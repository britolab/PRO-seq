#$ -S /bin/bash
#$ -N parse_vibrant_gbk
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/parse_vibrant_gbk_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/parse_vibrant_gbk_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=5G
#$ -q long.q@cbsubrito
#$ -t 1-2

# only for high quality integrated prophage

WRK=/workdir/users/acv46/stool_PROSeq2/phage
export OMP_NUM_THREADS=8

DESIGN_FILE=$WRK/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

base=$(echo ${NAME} | awk -F_ '{print $1}')

gbk=$WRK/vibrant/${base}/VIBRANT_${base}*/VIBRANT_phages_${base}*/*.phages_combined.gbk
pro=$WRK/vibrant/${base}/VIBRANT_${base}*/VIBRANT_results_${base}*/VIBRANT_integrated_prophage_coordinates_*.tsv

list=$(sed '1d' ${pro})

cd $WRK/vibrant

while read result; do

	contig=$(echo $result | awk '{print $1}')
	phage=$(echo $result | awk '{print $2}')
	start=$(echo $result | awk '{print $6}')
	end=$(echo $result | awk '{print $7}')

	genelist=$(grep "${phage}_" $gbk | sed 's|/locus_tag=||g' | tr -d '"' | tr -d [[:blank:]])

	while read gene; do

		gstart=$(grep -B1 "${gene}\"" $gbk | grep "CDS" | awk '{print $2}' | sed "s/complement(//g" | awk -F. '{print $1}')
		gend=$(grep -B1 "${gene}\"" $gbk | grep "CDS" | awk '{print $2}' | sed "s/complement(//g" | awk -F. '{print $3}' | sed "s/)//g")
		gname=$(grep -A2 "${gene}\"" $gbk | grep "/product" | cut -d \" -f2)
		test=$(grep -B1 "${gene}\"" $gbk | grep "CDS" | awk '{print $2}')

		if [ -z "$gname" ]; then
			gname=NA
		fi

		grep -B1 "${gene}\"" $gbk | grep "CDS" | awk '{print $2}' | grep -q "complement"
		if [ $? -eq 0 ]; then
			strand=-1
		else
			strand=1
		fi

		echo -e "$contig\t$phage\t$start\t$end\t$gene\t$gstart\t$gend\t$gname\t$strand" \
		>> ${base}_phage_summary.txt

	done <<< "$genelist"

done <<< "$list"

#grep -v "/gene" ${base}_temp.txt > ${base}_antismash_summary.txt
#rm ${base}_temp.txt
