#$ -S /bin/bash
#$ -N humann2
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/humann2_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/humann2_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=40G
#$ -q long.q@cbsubrito
#$ -t 1-10

WRK=/workdir/users/acv46/stool_PROSeq2
OUT=$WRK/out/humann2
MPH=$WRK/out/metaphlan
FQ=$WRK/fastq/clean
HUM=/home/britolab/acv46/.local/bin/humann2
export OMP_NUM_THREADS=8

# see $HUM -h for options

# download databases before running
# /home/britolab/acv46/.local/bin/humann2_databases \
# --download uniref uniref90_ec_filtered_diamond /workdir/refdbs/HUMAnN2/

LIST=$OUT/samples.txt
  DESIGN=$(sed -n "${SGE_TASK_ID}p" $LIST)
  NAME=`basename "$DESIGN"`

echo "${NAME} -- starting HUMAnN2"
echo "${NAME} -- checking for concatenated paired-end reads"

if [ ! -f $FQ/${NAME}.fastq.gz ]; then

	echo "${NAME} -- concatenating paired-end reads"

	cat $FQ/${NAME}_R1.fastq.gz $FQ/${NAME}_R2.fastq.gz \
	> $FQ/${NAME}.fastq.gz

fi

echo "${NAME} -- running HUMAnN2 with ${OMP_NUM_THREADS} threads"

# must use a version of diamond compatible with supplied databases

# resolve this
# https://groups.google.com/g/humann-users/c/MSVmXC7DSW0/m/jXl7VGZpCQAJ

$HUM \
	-i $FQ/${NAME}.fastq.gz \
	-o $OUT/${NAME} \
	--taxonomic-profile $MPH/${NAME}_profile.txt \
	--input-format fastq.gz \
	--threads ${OMP_NUM_THREADS} \
	--protein-database /workdir/refdbs/HUMAnN2/uniref \
	--nucleotide-database /workdir/refdbs/HUMAnN2/chocophlan \
	--metaphlan /programs/MetaPhlAn-2.0 \
	--diamond /programs/diamond-0.8.34

# diamond troubleshooting
# run /home/acv46/miniconda3/bin/diamond dbinfo --db <dmnd database>
# refer to https://github.com/bbuchfink/diamond/wiki/5.-Advanced-topics#database-format-versions

DB=uniref50
echo "${NAME} -- annotating output table with feature names from ${DB}"

${HUM}_rename_table \
	--input $OUT/${NAME}/${NAME}_genefamilies.tsv \
	--output $OUT/${NAME}/${NAME}_genefamilies-names.tsv \
	--names ${DB}

echo "${NAME} -- normalizing RPKs to relative abundance"

${HUM}_renorm_table \
	--input $OUT/${NAME}/${NAME}_genefamilies-names.tsv \
	--output $OUT/${NAME}/${NAME}_genefamilies-names-cpm.tsv \
	--units cpm \
	--update-snames

echo "${NAME} -- HUMAnN2 done"

echo "${NAME} -- removing concatenated reads to save space"

rm $FQ/${NAME}.fastq.gz
