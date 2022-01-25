#$ -S /bin/bash
#$ -N phaster
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq2/log/phaster_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq2/log/phaster_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq2/phage/phaster
#$ -l h_vmem=5G
#$ -q long.q@cbsubrito
#$ -t 1-2

# script to submit contigs to phaster API
# first, trim sequences shorter than 2000 nt
# second, split up fastas into files smaller than 25 Mb (max: 26)
# third, submit fasta pieces to phaster API

BBMAP=/programs/bbmap-38.86/reformat.sh
SQK=/programs/seqkit-0.15.0/seqkit
WRK=/workdir/users/acv46/stool_PROSeq2/phage/phaster
export OMP_NUM_THREADS=2

DESIGN_FILE=$WRK/fastas.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

cd $WRK

numseq1=$(grep ">" ${NAME}.fasta | wc -l)
echo "${NAME} -- $numseq1 sequences in input"

echo "${NAME} -- removing seqs shorter than 2kb"

$BBMAP \
	in=${NAME}.fasta \
	out=${NAME}_2kb.fasta \
	minlength=2000

numseq2=$(grep ">" ${NAME}_2kb.fasta | wc -l)
echo "${NAME} -- $numseq2 sequences after trimming"

size=$(stat -c%s ${NAME}.fasta)
max=25000000
numfiles=$(expr $size / $max + 1)

if (( size > max )); then

	echo "${NAME} -- fasta is larger than 25 Mb, splitting..."
	# split into $numfiles files
	$SQK split2 -l 20M ${NAME}.fasta
	# enter split dir and loop through files
	cd $WRK/${NAME}.fasta.split
	for file in *.fasta; do
		 wget --post-file="${file}" \
		"http://phaster.ca/phaster_api?contigs=1" \
		-O ${file}_phaster.out
		sleep 30
	done


## need to finish this
else

fi



