#$ -S /bin/bash
#$ -N pasrse_invertons
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/parse_invertons_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/parse_invertons_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=20G
#$ -q long.q@cbsubrito
#$ -t 1-4

# parse gtf (prokka) and PhaseFinder outputs to get invertase proximal to IRs
# search flanked sequence for promoters

INV=/workdir/users/acv46/stool_PROSeq2/invertons
PRK=/workdir/users/acv46/stool_PROSeq2/deseq/prokka
NNPP2=/home/acv46/MetaRon/NNPP2.2/bin/fa2TDNNpred.linux
gref=${INV}/ref_gtdbtk.txt
dref=${INV}/ref_dastool.txt

DESIGN_FILE=$INV/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

GTF=${PRK}/${NAME}/${NAME}_CUSTOM_FType.gtf
base=$(echo "${NAME}" | awk -F_ '{print $1}')


# get invertase genes from prokka gtf
grep -i "invertase" $GTF | \
	awk -F$'\t' '{print $1,$4,$5,$7,$9}' | \
	awk -F$';' '{print $1,$3}' | \
	sed 's/ \+/ /g' | tr [[:blank:]] '\t' | \
	awk -F$'\t' '{print $1,$2,$3,$4,$6,$8}' | \
	sed 's/\"//g' > ${INV}/${NAME}_invertase.txt



# loop through list of invertase genes in contigs
while read line; do

	# get inverted repeats on contig
	node=$(echo $line | awk '{print $1}')
	grep $node ${INV}/out_mismatch3/${base}_*.einverted.tab \
	> ${INV}/${node}_tmp

	# get bin ID and gtbdtk taxa
	ddir=$(grep "${base}" $dref)
	gdir=$(grep "${base}" $gref)
	bin=$(grep -r "${node}" $ddir | awk -F$'/' '{print $9}' | sed "s/.fa:.*//g")
	lineage=$(grep -w "^${bin}" $gdir | awk '{print $2}')

	# loop through inverted repeat hits
	while read tmp; do

		# create fasta for sequence flanked by IRs
		start=$(echo $tmp | awk '{print $3 + 1}')
		stop=$(echo $tmp | awk '{print $4 - 1}')
		echo $tmp | awk '{print $7}' | \
		sed "s/^/>${node}_${start}-${stop}\n/g" \
		> ${INV}/${node}_${start}-${stop}.fa

		# check fasta for putative promoters
		$NNPP2 -r -t 0.99 ${INV}/${node}_${start}-${stop}.fa \
		> ${INV}/${node}_${start}-${stop}.nnpp2

		# get the coordinates for the best promoter hit
		# append other info
		# write out
		best=$(grep -w "^Prediction:" ${INV}/${node}_${start}-${stop}.nnpp2 | sort -r | head -1)
		pseq=$(grep -B1 "${best}" ${INV}/${node}_${start}-${stop}.nnpp2 | head -1 | sed "s/.*: //g")
		grep -B2 "${best}" ${INV}/${node}_${start}-${stop}.nnpp2 | \
		head -1 | awk -F' ' '{print $4,$6}' | \
		sed 's/,//g' | sed "s/^/${bin}\t${lineage}\t${line}\t${tmp}\t${start}\t${stop}\t${pseq}\t/g" | \
		tr [[:blank:]] '\t' >> ${INV}/${NAME}.parsed.txt

		rm ${INV}/${node}_${start}-${stop}.fa
		rm ${INV}/${node}_${start}-${stop}.nnpp2


	done < ${INV}/${node}_tmp
	rm ${INV}/${node}_tmp

done < ${INV}/${NAME}_invertase.txt
