#$ -S /bin/bash
#$ -N pileup_du
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/pileup_du_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/pileup_du_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=100G
#$ -q long.q@cbsubrito
#$ -t 1-2

# see /workdir/users/acv46/EC_PROSeq/pileup/_README
# uses existing merged bams from get_pileup.sh
# subsets bams for MAPQ >= qscore and recounts pileups
# merges pileups into single file

export OMP_NUM_THREADS=20
BAM=/workdir/users/acv46/stool_PROSeq2/deseq/bam
OUT=$BAM/pileup

qscore=30

List=$BAM/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $List)
        NAME=`basename "$DESIGN"`

cd $OUT

#echo "${NAME} -- filtering for MAPQ >= ${qscore}"

# filter for alignments with following characteristics
## mapq score >= 30
## read paired and in proper pair
## read not unmapped, mate not unmapped, not not primary alignment, not supplementary alignment

#samtools view -@ $OMP_NUM_THREADS \
#	-f 3 -F 2316 \
#	-b -q ${qscore} \
#	-o ${NAME}_PRO_merge_q${qscore}.bam \
#	${NAME}_PRO_merge.bam

#samtools view -@ $OMP_NUM_THREADS \
#        -f 3 -F 2316 \
#	-b -q ${qscore} \
#       -o ${NAME}_RNA_merge_q${qscore}.bam \
#        ${NAME}_RNA_merge.bam

# USE EXISTING BAM FILES

echo "${NAME} -- generating stranded paired-end coverage for PRO-seq alignments"
echo "${NAME} --  ... full interval"

# no pc, no du, just visualize mapping
# might have to flip 5' and 3'
# swap plus and minus strand for corrent orientation
# 5' is actually 3' of pair?
# https://www.nature.com/articles/nprot.2016.086

# plus is minus, minus is plus
# 3 is 5, 5 is 3

bedtools genomecov -d -strand + -ibam ${NAME}_PRO_merge_q${qscore}.bam \
        > ${NAME}_PRO_minus_q${qscore}_full.cov
bedtools genomecov -d -strand - -ibam ${NAME}_PRO_merge_q${qscore}.bam \
        > ${NAME}_PRO_plus_q${qscore}_full.cov

echo "${NAME} --  ... 3' ends"

bedtools genomecov -5 -d -strand + -ibam ${NAME}_PRO_merge_q${qscore}.bam \
        > ${NAME}_PRO_minus_q${qscore}_3p.cov
bedtools genomecov -5 -d -strand - -ibam ${NAME}_PRO_merge_q${qscore}.bam \
        > ${NAME}_PRO_plus_q${qscore}_3p.cov

echo "${NAME} --  ... 5' ends"

bedtools genomecov -3 -d -strand + -ibam ${NAME}_PRO_merge_q${qscore}.bam \
	> ${NAME}_PRO_minus_q${qscore}_5p.cov
bedtools genomecov -3 -d -strand - -ibam ${NAME}_PRO_merge_q${qscore}.bam \
	> ${NAME}_PRO_plus_q${qscore}_5p.cov

echo "${NAME} -- generating stranded paired-end coverage for RNA-seq alignments"
echo "${NAME} --  ... full interval"

bedtools genomecov -d -strand + -ibam ${NAME}_RNA_merge_q${qscore}.bam \
        > ${NAME}_RNA_minus_q${qscore}_full.cov
bedtools genomecov -d -strand - -ibam ${NAME}_RNA_merge_q${qscore}.bam \
        > ${NAME}_RNA_plus_q${qscore}_full.cov

echo "${NAME} --  ... 5' ends"

bedtools genomecov -5 -d -strand + -ibam ${NAME}_RNA_merge_q${qscore}.bam \
        > ${NAME}_RNA_minus_q${qscore}_3p.cov
bedtools genomecov -5 -d -strand - -ibam ${NAME}_RNA_merge_q${qscore}.bam \
        > ${NAME}_RNA_plus_q${qscore}_3p.cov

echo "${NAME} --  ... 3' ends"

bedtools genomecov -du -3 -d -strand + -ibam ${NAME}_RNA_merge_q${qscore}.bam \
        > ${NAME}_RNA_minus_q${qscore}_5p.cov
bedtools genomecov -du -3 -d -strand - -ibam ${NAME}_RNA_merge_q${qscore}.bam \
        > ${NAME}_RNA_plus_q${qscore}_5p.cov

for file in ${NAME}_*_q${qscore}_*.cov; do

	name=$(echo "$file" | sed "s/.cov//g")
	sed -i 1i"contig\tposition\t$name" $file

done

paste ${NAME}_*_q${qscore}_*.cov | \
	awk '{print $1,$2,$3,$6,$9,$12,$15,$18,$21,$24,$27,$30,$33,$36}' | \
	tr [[:blank:]] "\t" \
	> ${NAME}_coverage_q${qscore}_flipEnd_flipStrand.txt

echo "${NAME} -- done!"

# remove intermediate files
rm paste ${NAME}_*_q${qscore}_*.cov
#rm ${NAME}_PRO_merge_q${qscore}.bam
#rm ${NAME}_RNA_merge_q${qscore}.bam
