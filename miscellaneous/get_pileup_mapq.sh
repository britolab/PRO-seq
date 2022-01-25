#$ -S /bin/bash
#$ -N pileup_mapq
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/pileup20_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/pileup20_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=100G
#$ -q long.q@cbsubrito
#$ -t 1-2

# see /workdir/users/acv46/EC_PROSeq/pileup/_README
# uses existing merged bams from get_pileup.sh
# subsets bams for MAPQ >= qscore and recounts pileups
# merges pileups into single file

export OMP_NUM_THREADS=16
BAM=/workdir/users/acv46/stool_PROSeq2/deseq/bam
OUT=$BAM/pileup

qscore=30

List=$BAM/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $List)
        NAME=`basename "$DESIGN"`

cd $OUT

echo "${NAME} -- filtering for MAPQ >= ${qscore}"

# filter for alignments with following characteristics
## mapq score >= 30
## read paired and in proper pair
## read not unmapped, mate not unmapped, not not primary alignment, not supplementary alignment

samtools view -@ $OMP_NUM_THREADS \
	-f 3 -F 2316 \
	-b -q ${qscore} \
	-o ${NAME}_PRO_merge_q${qscore}.bam \
	${NAME}_PRO_merge.bam

samtools view -@ $OMP_NUM_THREADS \
        -f 3 -F 2316 \
	-b -q ${qscore} \
        -o ${NAME}_RNA_merge_q${qscore}.bam \
        ${NAME}_RNA_merge.bam

echo "${NAME} -- generating stranded paired-end coverage for PRO-seq alignments, 3' only"

bedtools genomecov -pc -3 -d -strand + -ibam ${NAME}_PRO_merge_q${qscore}.bam \
	> ${NAME}_PRO_plus_q${qscore}_pp.cov
bedtools genomecov -pc -3 -d -strand - -ibam ${NAME}_PRO_merge_q${qscore}.bam \
	> ${NAME}_PRO_minus_q${qscore}_pp.cov

echo "${NAME} -- generating paired-end coverage for RNA-seq alignments"

bedtools genomecov -pc -d -strand + -ibam ${NAME}_RNA_merge_q${qscore}.bam \
	> ${NAME}_RNA_plus_q${qscore}_pp.cov
bedtools genomecov -pc -d -strand - -ibam ${NAME}_RNA_merge_q${qscore}.bam \
	> ${NAME}_RNA_minus_q${qscore}_pp.cov

for file in ${NAME}_*_q${qscore}_pp.cov; do

	name=$(echo "$file" | sed "s/.cov//g")
	sed -i 1i"contig\tposition\t$name" $file

done

paste ${NAME}_*_q${qscore}_pp.cov | \
	awk '{print $1,$2,$3,$6,$9,$12}' | \
	tr [[:blank:]] "\t" \
	> ${NAME}_coverage_q${qscore}_pp.txt

echo "${NAME} -- done!"

# remove intermediate files
rm paste ${NAME}_*_q${qscore}_pp.cov
rm ${NAME}_PRO_merge_q${qscore}.bam
rm ${NAME}_RNA_merge_q${qscore}.bam
