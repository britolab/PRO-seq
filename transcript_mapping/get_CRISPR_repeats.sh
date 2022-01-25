#$ -S /bin/bash
#$ -N CRISPR_repeats
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq2/log/CRISPR_repeats_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq2/log/CRISPR_repeats_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=50G
#$ -q long.q@cbsubrito
#$ -t 1-4

# takes minced output from prokka run and blasts against contigs

WRK=/workdir/users/acv46/stool_PROSeq2
BLASTDB=$WRK/phage/virsorter
OUT=$WRK/phage/CRISPR
CONTIG=$WRK/deseq/ref
GTF=$WRK/deseq/prokka

samtools=/programs/samtools-1.11/bin/samtools
seqtk=/programs/seqtk/seqtk
blastn=/programs/ncbi-blast-2.9.0+/bin/blastn
export OMP_NUM_THREADS=4

DESIGN_FILE=$OUT/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

mkdir -p $OUT/${NAME}
cd $OUT/${NAME}

base=$(echo $NAME | awk -F_ '{print $1}')
db=$(grep "${base}" ${OUT}/dblist.txt)

awk '$3 == "repeat_region" {print $0}' $GTF/${NAME}/${NAME}_CUSTOM_FType.gtf \
	> minced_repeat_regions.gtf

while read line; do

	node=$(echo "${line}" | awk -F$'\t' '{print $1}')
	begin=$(echo "${line}" | awk -F$'\t' '{print $4}')
	end=$(echo "${line}" | awk -F$'\t' '{print $5}')
	repeat=$(echo "$line" | awk -F$'\t' '{print $9}' | awk -F'"' '{print $8}')

	$samtools faidx $CONTIG/${NAME}.short.fasta ${node}:${begin}-${end} | \
		$seqtk seq - | \
		sed "s/${repeat}/|/g" | \
		tr '|' '\n' | sed "/^$/d; 1d; h; s/.*//; s/$/>${node}_CRISPR_${begin}-${end}/; G" | \
		awk '{if ($0 ~/^>/) {h[$1]++; $1=$1 "_spacer" h[$1]} print}' | \
		$blastn -query - -db ${BLASTDB}/${db}/blastdb/${base}_phage -outfmt 6 \
		>> ${NAME}_CRISPR_blastn.txt

done < minced_repeat_regions.gtf
