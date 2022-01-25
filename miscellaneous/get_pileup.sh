#$ -S /bin/bash
#$ -N pileup
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/pileup_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/pileup_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=100G
#$ -q long.q@cbsubrito
#$ -t 1-2

# see /workdir/users/acv46/EC_PROSeq/pileup/_README

export OMP_NUM_THREADS=16
BAM=/workdir/users/acv46/stool_PROSeq2/deseq/bam
OUT=$BAM/pileup

List=$BAM/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $List)
        NAME=`basename "$DESIGN"`

cd $BAM

samtools merge -@ $OMP_NUM_THREADS \
	$OUT/${NAME}_PRO_merge.bam \
	${NAME}_PRO1.allcontigs.bam \
	${NAME}_PRO2.allcontigs.bam

samtools merge -@ $OMP_NUM_THREADS \
        $OUT/${NAME}_RNA_merge.bam \
        ${NAME}_RNA1.allcontigs.bam \
        ${NAME}_RNA2.allcontigs.bam

cd $OUT

bedtools genomecov -3 -d -strand + -ibam ${NAME}_PRO_merge.bam > ${NAME}_PRO_plus.cov
bedtools genomecov -3 -d -strand - -ibam ${NAME}_PRO_merge.bam > ${NAME}_PRO_minus.cov
bedtools genomecov -d -strand + -ibam ${NAME}_RNA_merge.bam > ${NAME}_RNA_plus.cov
bedtools genomecov -d -strand - -ibam ${NAME}_RNA_merge.bam > ${NAME}_RNA_minus.cov

#for file in *.cov; do name=$(echo "$file" | sed "s/.cov//g"); sed -i 1i"genome\tposition\t$name" $file ; done

#paste ./*.cov | awk '{print $1,$2,$3,$6,$9,$12,$15,$18,$21,$24}' | tr [[:blank:]] "\t" > all_coverage.txt
