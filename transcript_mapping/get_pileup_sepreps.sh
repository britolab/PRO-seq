#$ -S /bin/bash
#$ -N pileup
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/pileup_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/pileup_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=40G
#$ -q long.q@cbsubrito
#$ -t 1-8

# see /workdir/users/acv46/EC_PROSeq/pileup/_README
# redo to determine which files are mapping to wrong strand?
# seperate replicates, not merged

export OMP_NUM_THREADS=8
BAM=/workdir/users/acv46/stool_PROSeq2/deseq/bam
OUT=$BAM/pileup
qscore=30

List=$BAM/samples_reps.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $List)
        NAME=`basename "$DESIGN"`

cd $BAM

echo "${NAME} -- filtering for MAPQ >= ${qscore}"

samtools view -@ $OMP_NUM_THREADS \
        -b -q ${qscore} \
        -o $OUT/${NAME}_q${qscore}.bam \
        ${NAME}.bam

cd $OUT

if [[ ${NAME} == *"PRO"* ]]; then

	echo "${NAME} -- PRO-seq sample, generating stranded 3' coverage files"

	bedtools genomecov -3 -d -strand + -ibam ${NAME}_q${qscore}.bam \
		> ${NAME}_q${qscore}_plus.cov
	bedtools genomecov -3 -d -strand - -ibam ${NAME}_q${qscore}.bam \
                > ${NAME}_q${qscore}_minus.cov

elif [[ ${NAME} == *"RNA"* ]]; then

	echo "${NAME} -- RNA-seq sample, generating stranded coverage files"

	bedtools genomecov -d -strand + -ibam ${NAME}_q${qscore}.bam \
                > ${NAME}_q${qscore}_plus.cov
        bedtools genomecov -d -strand - -ibam ${NAME}_q${qscore}.bam \
                > ${NAME}_q${qscore}_minus.cov

fi

echo "${NAME} -- done!"

#for file in *.cov; do name=$(echo "$file" | sed "s/.cov//g"); sed -i 1i"genome\tposition\t$name" $file ; done

#paste ./*.cov | awk '{print $1,$2,$3,$6,$9,$12,$15,$18,$21,$24}' | tr [[:blank:]] "\t" > all_coverage.txt
