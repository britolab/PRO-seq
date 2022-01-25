#$ -S /bin/bash
#$ -N pileup_correct
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/pileup_correct_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/pileup_correct_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=50G
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

# start with pre-filtered bam file

# filter for alignments with following characteristics
## mapq score >= 30
## read paired and in proper pair
## read not unmapped, mate not unmapped, not not primary alignment, not supplementary alignment
# -f 3 -F 2316

# P5 Illumina adapter on 3' end of nascent transcript
# 5' end of first in pair gives 3' end of read
# 5' end of second in pair gives 5' end of read
# flip strands

# sort the files before combining, dumbass

echo "${NAME} -- generating stranded paired-end coverage for PRO-seq alignments"
echo "${NAME} --  ... full interval"

bedtools genomecov -pc -d -strand + -ibam ${NAME}_PRO_merge_q${qscore}.bam | \
	sort -nk1.6 -nk2 \
        > ${NAME}_PRO_minus_q${qscore}_full.cov
echo "${NAME}_PRO_minus_q${qscore}_full.cov -- $(wc -l ${NAME}_PRO_minus_q${qscore}_full.cov)"
bedtools genomecov -pc -d -strand - -ibam ${NAME}_PRO_merge_q${qscore}.bam | \
	sort -nk1.6 -nk2 \
        > ${NAME}_PRO_plus_q${qscore}_full.cov
echo "${NAME}_PRO_plus_q${qscore}_full.cov -- $(wc -l ${NAME}_PRO_plus_q${qscore}_full.cov)"

echo "${NAME} --  ... 3' fragment ends"

samtools view -@ $OMP_NUM_THREADS -f 64 -b ${NAME}_PRO_merge_q${qscore}.bam | \
	bedtools genomecov -5 -d -strand + -ibam - | \
        sort -nk1.6 -nk2 \
        > ${NAME}_PRO_minus_q${qscore}_3p.cov
echo "${NAME}_PRO_minus_q${qscore}_3p.cov -- $(wc -l ${NAME}_PRO_minus_q${qscore}_3p.cov)"
samtools view -@ $OMP_NUM_THREADS -f 64 -b ${NAME}_PRO_merge_q${qscore}.bam | \
        bedtools genomecov -5 -d -strand - -ibam - | \
        sort -nk1.6 -nk2 \
        > ${NAME}_PRO_plus_q${qscore}_3p.cov
echo "${NAME}_PRO_plus_q${qscore}_3p.cov -- $(wc -l ${NAME}_PRO_plus_q${qscore}_3p.cov)"

echo "${NAME} --  ... 5' fragment ends"

# because proper pairs map to opposite strands, 5' fragment end is correct strand

samtools view -@ $OMP_NUM_THREADS -f 128 -b ${NAME}_PRO_merge_q${qscore}.bam | \
        bedtools genomecov -5 -d -strand + -ibam - | \
        sort -nk1.6 -nk2 \
        > ${NAME}_PRO_plus_q${qscore}_5p.cov
echo "${NAME}_PRO_plus_q${qscore}_5p.cov -- $(wc -l ${NAME}_PRO_plus_q${qscore}_5p.cov)"
samtools view -@ $OMP_NUM_THREADS -f 128 -b ${NAME}_PRO_merge_q${qscore}.bam | \
        bedtools genomecov -5 -d -strand - -ibam - | \
        sort -nk1.6 -nk2 \
        > ${NAME}_PRO_minus_q${qscore}_5p.cov
echo "${NAME}_PRO_minus_q${qscore}_5p.cov -- $(wc -l ${NAME}_PRO_minus_q${qscore}_5p.cov)"

echo "${NAME} -- generating stranded paired-end coverage for RNA-seq alignments"
echo "${NAME} --  ... full interval"

bedtools genomecov -pc -d -strand + -ibam ${NAME}_RNA_merge_q${qscore}.bam | \
        sort -nk1.6 -nk2 \
        > ${NAME}_RNA_plus_q${qscore}_full.cov
echo "${NAME}_RNA_plus_q${qscore}_full.cov -- $(wc -l ${NAME}_RNA_plus_q${qscore}_full.cov)"
bedtools genomecov -pc -d -strand - -ibam ${NAME}_RNA_merge_q${qscore}.bam | \
        sort -nk1.6 -nk2 \
        > ${NAME}_RNA_minus_q${qscore}_full.cov
"${NAME}_RNA_minus_q${qscore}_full.cov -- $(wc -l ${NAME}_RNA_minus_q${qscore}_full.cov)"

for file in ${NAME}_*_q${qscore}_*.cov; do

	name=$(echo "$file" | sed "s/.cov//g")
	sed -i 1i"contig\tposition\t$name" $file

done

paste ${NAME}_*_q${qscore}_*.cov | \
	awk '{print $1,$2,$3,$6,$9,$12,$15,$18,$21,$24}' | \
	tr [[:blank:]] "\t" \
	> ${NAME}_coverage_q${qscore}_correct3.txt

echo "${NAME} -- done!"

# remove intermediate files
rm ${NAME}_*_q${qscore}_*.cov
#rm ${NAME}_PRO_merge_q${qscore}.bam
#rm ${NAME}_RNA_merge_q${qscore}.bam
