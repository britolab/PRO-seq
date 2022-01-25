#$ -S /bin/bash
#$ -N iRep
#$ -V
#$ -o /workdir/users/acv46/mgmAssembly/log/iRep_$JOB_ID.out
#$ -e /workdir/users/acv46/mgmAssembly/log/iRep_$JOB_ID.err
#$ -wd /workdir/users/acv46/mgmAssembly
#$ -l h_vmem=40G
#$ -t 1-3
#$ -q long.q@cbsubrito

#uses bam alignments and cleaned bins generated with getBins.sh

WRK=/workdir/users/acv46/CDC_Sep2019/iRep
BAM=$WRK/bam
BIN=$WRK/cleanedBins
IREP=/home/acv46/miniconda3/bin/iRep
OUT=$WRK/out
SHRINK=/home/britolab/acv46/shrinksam/shrinksam-master/shrinksam

#Full directory to list of sample names
list=$WRK/samples.txt

List=${list} #The location and name of your "design file" or list of file names in a tex$
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $List)
        samp=`basename "$DESIGN"`

echo "Converting ${samp}.sorted.bam to sam"
samtools view -h -o $BAM/${samp}.sorted.sam $BAM/${samp}.sorted.bam

echo "Removing sequences in ${samp}.sorted.sam that failed to map"
$SHRINK -i $BAM/${samp}.sorted.sam -k $BAM/${samp}.shrink.sam

mkdir -p $OUT/${samp}

awk '{print $1}' $BIN/${samp}.goodBins.summary.txt > $OUT/${samp}/${samp}.bins

countbins=$(wc -l $OUT/${samp}/${samp}.bins | awk -F' ' '{print $1}')
echo "There are ${countbins} HQ bins for ${samp}"

while read line; do

	donebins=$(grep -n "${line}" $OUT/${samp}/${samp}.bins | \
		awk -F: '{print $1}')

	echo "Processing ${samp} HQ bin ${line} (${donebins} of ${countbins})"
	$IREP -f $WRK/das/${samp}/${samp}_DASTool_bins/${line}.contigs.fa \
	-s $BAM/${samp}.shrink.sam \
	-o $OUT/${samp}/${samp}.${line}.iRep \
	-t 12

	kraken=$(cat $WRK/das/${samp}/kraken/${line}.contigs.fa.kraken.weighted_besttaxid.txt | \
		awk -F$'\t' '{print $8}')
	echo "kraken annotation for ${samp} bin ${line}: ${kraken}"

	checkm=$(grep "${line}.contigs" $WRK/das/${samp}/checkm_lineage/${samp}.stats | \
		awk -F$'\t' '{print $2}')
	echo "checkm annotation for ${samp} bin ${line}: ${checkm}"

	echo -e "${samp}\t${line}\t${kraken}\t${checkm}" >> $OUT/summary.$JOB_ID

	echo "Finished ${samp} HQ bin ${line}"

done < $OUT/${samp}/${samp}.bins

echo "${samp} DONE"
echo "removing sam files"

rm $BAM/${samp}.sorted.sam
rm $BAM/${samp}.shrink.sam
rm $OUT/${samp}/${samp}.bins
