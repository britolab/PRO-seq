#$ -S /usr/bin/env Rscript
#$ -N filter_bam
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/filter_bam_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/filter_bam_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=60G
#$ -q long.q@cbsubrito
#$ -t 1-4

# see /workdir/users/acv46/EC_PROSeq/pileup/_README
# uses existing merged bams from get_pileup.sh
# subsets bams for MAPQ >= qscore and recounts pileups
# merges pileups into single file

export OMP_NUM_THREADS=8
BAM=/workdir/users/acv46/stool_PROSeq2/deseq/bam
OUT=$BAM/pileup

qscore=30

List=$OUT/bams.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $List)
        NAME=`basename "$DESIGN"`

cd $OUT

base=$(echo ${NAME} | awk -F. '{print $1}')
echo "${base} -- filtering for MAPQ >= ${qscore}"

# filter for alignments with following characteristics
## mapq score >= 30
## read paired and in proper pair
## read not unmapped, mate not unmapped, not not primary alignment, not supplementary alignment

samtools view -@ $OMP_NUM_THREADS \
	-f 3 -F 2316 \
	-b -q ${qscore} \
	-o ${base}_q${qscore}.bam \
	${NAME}

echo "${base} -- done"
