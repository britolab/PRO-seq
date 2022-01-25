#$ -S /bin/bash
#$ -N featureCounts
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq2/log/featureCounts_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq2/log/featureCounts_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=60G
#$ -q long.q@cbsubrito
#$ -t 1-4

# iterates over different GTF-bam pairings
# runs REV mapping (-s 2) for both RNA-seq and PRO-seq

WRK=/workdir/users/acv46/stool_PROSeq2/deseq
BAM=$WRK/bam
FCT=/home/britolab/acv46/subread-2.0.2-source/bin/featureCounts
export OMP_NUM_THREADS=8

DESIGN_FILE=$WRK/jobs.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

samp=$(echo $NAME | awk -F_ '{print $1}')
lib=$(echo $NAME | awk -F_ '{print $2}')

OUT=$WRK/featureCounts/${NAME}
mkdir -p $OUT

# paired reads (-p)
# only count read pairs where both ends are aligned (-B)
# do not count read pairs if reads map to different contigs or different strands (-C)
# perform read counting at feature level (-f)

echo "${NAME} -- filtering bams by MAPQ, etc."

samtools view -@ $OMP_NUM_THREADS \
        -f 3 -F 2316 \
        -b -q 30 \
        -o $BAM/${samp}_PRO1.${lib}_q30.bam \
        $BAM/${samp}_PRO1.${lib}.bam

samtools view -@ $OMP_NUM_THREADS \
        -f 3 -F 2316 \
        -b -q 30 \
        -o $BAM/${samp}_PRO2.${lib}_q30.bam \
        $BAM/${samp}_PRO2.${lib}.bam

samtools view -@ $OMP_NUM_THREADS \
        -f 3 -F 2316 \
        -b -q 30 \
        -o $BAM/${samp}_RNA1.${lib}_q30.bam \
        $BAM/${samp}_RNA1.${lib}.bam

samtools view -@ $OMP_NUM_THREADS \
        -f 3 -F 2316 \
        -b -q 30 \
        -o $BAM/${samp}_RNA2.${lib}_q30.bam \
        $BAM/${samp}_RNA2.${lib}.bam

# get sense (reverse) counts

echo "${NAME} -- getting sense feature counts"

$FCT \
	-a $WRK/prokka/${NAME}/${NAME}_CUSTOM_FType.gtf \
	-o $OUT/${NAME}_SENSE_q30.txt \
	-F GTF \
	-t CDS,misc_RNA,repeat_region,rRNA,tmRNA,tRNA \
	-g feature_type \
	--extraAttributes Name,product,gene,Note,ID,Parent \
	-s 2 \
	--countReadPairs \
	-f -p -B -C \
	-T $OMP_NUM_THREADS \
	--verbose \
	$BAM/${samp}_PRO1.${lib}_q30.bam \
	$BAM/${samp}_PRO2.${lib}_q30.bam \
	$BAM/${samp}_RNA1.${lib}_q30.bam \
	$BAM/${samp}_RNA2.${lib}_q30.bam

# get antisense (forward) counts

echo "${NAME} -- getting antisense feature counts"

$FCT \
        -a $WRK/prokka/${NAME}/${NAME}_CUSTOM_FType.gtf \
        -o $OUT/${NAME}_ANTISENSE_q30.txt \
        -F GTF \
        -t CDS,misc_RNA,repeat_region,rRNA,tmRNA,tRNA \
        -g feature_type \
        --extraAttributes Name,product,gene,Note,ID,Parent \
        -s 1 \
        --countReadPairs \
        -f -p -B -C \
        -T $OMP_NUM_THREADS \
        --verbose \
        $BAM/${samp}_PRO1.${lib}_q30.bam \
        $BAM/${samp}_PRO2.${lib}_q30.bam \
        $BAM/${samp}_RNA1.${lib}_q30.bam \
        $BAM/${samp}_RNA2.${lib}_q30.bam

echo "${NAME} -- removing intermediate bams"

rm $BAM/${samp}_*.${lib}_q30.bam
