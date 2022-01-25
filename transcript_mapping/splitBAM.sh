#$ -S /bin/bash
#$ -N splitBAM
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/splitBAM_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/splitBAM_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=30G
#$ -q long.q@cbsubrito
#$ -t 1-16

# Splits BAM files into reads mapping to forward and reverse strands
# assumes PRO and RNA seq reads are reversed, so reads identified mapped to "reverse" strand are in fact forward
# http://seqanswers.com/forums/showthread.php?t=29399

WRK=/workdir/users/acv46/stool_PROSeq2/deseq/bam
DEP=$WRK/depth
mkdir -p $DEP

export OMP_NUM_THREADS=8
export PATH=/programs/samtools-1.11/bin:$PATH
export PATH=/programs/bedtools2-2.29.2/bin:$PATH


LIST=$WRK/bams.txt
  DESIGN=$(sed -n "${SGE_TASK_ID}p" $LIST)
  NAME=`basename "$DESIGN"`

BAM=$WRK/${NAME}.bam

# explain flags
## -f 147 -- read paired, mapped in proper pair, on reverse strand, second in pair
## -f 99 -- read paired, mapped in proper pair, mate on reverse strand, first in pair
## -f 83 -- read paired, mapped in proper pair, on reverse strand, first in pair
## -f 163 -- read paired, mapped in proper pair, mate on reverse strand, second in pair
## -F 2828 -- read unmapped, mate unmapped, not primary alignment, supplmentary alignment, read fails quality checks

# this script assumes input bams are sorted by position
# bedtools will not work otherwise

if (echo $NAME | grep -q "PRO"); then

	echo "processing ${NAME} -- 3' ends only"

	samtools view -b -@ $OMP_NUM_THREADS -f 147 -F 2828 $BAM > $DEP/${NAME}_147.tmp
	samtools view -b -@ $OMP_NUM_THREADS -f 99 -F 2828 $BAM > $DEP/${NAME}_99.tmp
	samtools view -b -@ $OMP_NUM_THREADS -f 83 -F 2828 $BAM > $DEP/${NAME}_83.tmp
	samtools view -b -@ $OMP_NUM_THREADS -f 163 -F 2828 $BAM > $DEP/${NAME}_163.tmp

	samtools merge -@ $OMP_NUM_THREADS $DEP/${NAME}_fwd.bam $DEP/${NAME}_147.tmp $DEP/${NAME}_99.tmp
	rm $DEP/${NAME}_147.tmp
	rm $DEP/${NAME}_99.tmp

	samtools merge -@ $OMP_NUM_THREADS $DEP/${NAME}_rev.bam $DEP/${NAME}_83.tmp $DEP/${NAME}_163.tmp
	rm $DEP/${NAME}_83.tmp
	rm $DEP/${NAME}_163.tmp

	bedtools genomecov -d -3 -ibam $DEP/${NAME}_fwd.bam > $DEP/${NAME}_fwd.txt
	echo "created ${NAME}_fwd.txt for 3' end mapping on + strand"

	bedtools genomecov -d -3 -ibam $DEP/${NAME}_rev.bam > $DEP/${NAME}_rev.txt
	echo "created ${NAME}_rev.txt for 3' end mapping on - strand"

	rm $DEP/${NAME}_fwd.bam
	rm $DEP/${NAME}_rev.bam
	echo "${NAME} -- removed intermediate bam files"

else

	echo "processing ${NAME} -- full interval"

        samtools view -b -@ $OMP_NUM_THREADS -f 147 -F 2828 $BAM > $DEP/${NAME}_147.tmp
        samtools view -b -@ $OMP_NUM_THREADS -f 99 -F 2828 $BAM > $DEP/${NAME}_99.tmp
        samtools view -b -@ $OMP_NUM_THREADS -f 83 -F 2828 $BAM > $DEP/${NAME}_83.tmp
        samtools view -b -@ $OMP_NUM_THREADS -f 163 -F 2828 $BAM > $DEP/${NAME}_163.tmp

        samtools merge -@ $OMP_NUM_THREADS $DEP/${NAME}_fwd.bam $DEP/${NAME}_147.tmp $DEP/${NAME}_99.tmp
        rm $DEP/${NAME}_147.tmp
        rm $DEP/${NAME}_99.tmp

        samtools merge -@ $OMP_NUM_THREADS $DEP/${NAME}_rev.bam $DEP/${NAME}_83.tmp $DEP/${NAME}_163.tmp
        rm $DEP/${NAME}_83.tmp
        rm $DEP/${NAME}_163.tmp

        bedtools genomecov -d -ibam $DEP/${NAME}_fwd.bam > $DEP/${NAME}_fwd.txt
        echo "created ${NAME}_fwd.txt for full interval mapping on + strand"

        bedtools genomecov -d -ibam $DEP/${NAME}_rev.bam > $DEP/${NAME}_rev.txt
        echo "created ${NAME}_rev.txt for full interval mapping on - strand"

        rm $DEP/${NAME}_fwd.bam
        rm $DEP/${NAME}_rev.bam
        echo "${NAME} -- removed intermediate bam files"

fi

#paste $DEP/${NAME}_fwd.txt <(awk -v add=$(awk '{ sum += $3 } END { print sum }' $DEP/${NAME}_fwd.txt) -F$'\t' '{print ($3 / add)*100000000}' $DEP/${NAME}_fwd.txt) | \
#awk -F$'\t' '{print $2,$4}' | \
#sed "1i position\t${NAME}" | \
#tr [[:blank:]] '\t' > $DEP/${NAME}_fwd_norm.txt

#paste $DEP/${NAME}_rev.txt <(awk -v add=$(awk '{ sum += $3 } END { print sum }' $DEP/${NAME}_rev.txt) -F$'\t' '{print ($3 / add)*100000000}' $DEP/${NAME}_rev.txt) | \
#awk -F$'\t' '{print $2,$4}' | \
#sed "1i position\t${NAME}" | \
#tr [[:blank:]] '\t' > $DEP/${NAME}_rev_norm.txt
