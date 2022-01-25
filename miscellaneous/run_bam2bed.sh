#$ -S /bin/bash
#$ -N bam2bed
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq2/log/bam2bed_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq2/log/bam2bed_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=10G
#$ -q long.q@cbsubrito
#$ -t 1-16

# convert bam alignments to plus- and minus-strand bed files
# report depth of 3' positions as bedgraph files

WRK=/workdir/users/acv46/stool_PROSeq2/deseq
BAM=$WRK/bam
BED=$WRK/bed
export OMP_NUM_THREADS=6
export PATH=/programs/bedtools2-2.29.2/bin:$PATH

DESIGN_FILE=$BAM/bams.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

CHRM=$(echo $NAME | sed 's/_PRO1\|_PRO2\|_RNA1\|_RNA2//g' | sed 's/\./_/g')

bamToBed -i $BAM/${NAME}.bam | \
	grep -P '\t\+$' | \
	sort -k 1,1 > $BED/${NAME}.plus.bed
bamToBed -i $BAM/${NAME}.bam | \
	grep -P '\t-$' | \
	sort -k 1,1 > $BED/${NAME}.minus.bed

# if PRO-seq, just get 3' end coverage
# if RNA-seq, get coverage over entire interval

if (echo $NAME | grep -q "PRO"); then

	echo "processing ${NAME} -- 3' ends only"

	genomeCoverageBed -i $BED/${NAME}.plus.bed \
		-g $WRK/ref/${CHRM}.short.fasta.chrmSize \
		-3 -bg > $BED/${NAME}.plus.bedgraph

	genomeCoverageBed -i $BED/${NAME}.minus.bed \
        	-g $WRK/ref/${CHRM}.short.fasta.chrmSize \
        	-3 -bg > $BED/${NAME}.minus.bedgraph

else

	echo "processing ${NAME} -- full interval"

        genomeCoverageBed -i $BED/${NAME}.plus.bed \
                -g $WRK/ref/${CHRM}.short.fasta.chrmSize \
                -bg > $BED/${NAME}.plus.bedgraph

        genomeCoverageBed -i $BED/${NAME}.minus.bed \
                -g $WRK/ref/${CHRM}.short.fasta.chrmSize \
                -bg > $BED/${NAME}.minus.bedgraph

fi


