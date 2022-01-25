#$ -S /bin/bash
#$ -N annotateCRISPR
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq2/log/annotateCRISPR_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq2/log/annotateCRISPR_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=10G
#$ -q long.q@cbsubrito
#$ -t 1-4

# works on blast output of get_CRISPR_repeats.sh
## get bin associated with contig (node)
## get checkm annotation / gtdbtk annotation of bin
## get coverage of contig
## get PRO-seq and RNA-seq expression across CRISPR array
## get bin associated with phage node (if it exists)
## get coverage of phage node
## get expression across phage node

WRK=/workdir/users/acv46/stool_PROSeq2
BLASTDB=$WRK/phage/virsorter
OUT=$WRK/phage/CRISPR
CONTIG=$WRK/deseq/ref
GTF=$WRK/deseq/prokka

DESIGN_FILE=$OUT/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

blastout=$OUT/${NAME}/${NAME}_CRISPR_blastn.txt
cd $OUT/${NAME}

base=$(echo $NAME | awk -F_ '{print $1}')
db=$(grep "${base}" ${OUT}/dblist2.txt)

while read line; do

	checkm_cov=/workdir/users/acv46/mgmAssembly/${db}_CAB/DASTool/checkM/${db}_checkM_coverage.tsv
	checkm_lin=/workdir/users/acv46/mgmAssembly/${db}_CAB/DASTool/checkM/${db}_checkM_lineage.txt
	fcounts=/workdir/users/acv46/stool_PROSeq2/deseq/featureCounts/${NAME}/${NAME}_SENSE.txt

	hostnode=$(echo "${line}" | awk -F$'\t' '{print $1}' | awk -F$'_CRISPR' '{print $1}')
	hostnodelen=$(echo "${line}" | awk -F_ '{print $4}')
	cstart=$(echo "${line}" | awk -F_ '{print $6}' | awk -F- '{print $1}')
	cend=$(echo "${line}" | awk -F_ '{print $6}' | awk -F- '{print $2}')
	phagenode=$(echo "${line}" | awk -F$'\t' '{print $2}' | awk -F$'_cov' '{print $1}')
	phagenodelen=$(echo "${line}" | awk -F$'\t' '{print $2}' | awk -F_ '{print $4}')
	pstart=$(echo "${line}" | awk -F$'\t' '{print $9}')
	pend=$(echo "${line}" | awk -F$'\t' '{print $10}')
	blastscore=$(echo "${line}" | awk -F$'\t' '{print $3}')
	hostbin=$(grep "${hostnode}" $checkm_cov | awk -F$'\t' '{print $2}')
	hostannot=$(grep -w "${hostbin}" $checkm_lin | awk '{print $2}')
	hostcov=$(grep "${hostnode}" $checkm_cov | awk -F$'\t' '{print $5}')
	phagebin=$(grep "${phagenode}" $checkm_cov | awk -F$'\t' '{print $2}')
	phageannot=$(grep -w "${phagebin}" $checkm_lin | awk '{print $2}')
	phagecov=$(grep "${phagenode}" $checkm_cov | awk -F$'\t' '{print $5}')
	hostpro1=$(grep -P "repeat_region\t${hostnode}" $fcounts | awk '{print $13}')
	hostpro2=$(grep -P "repeat_region\t${hostnode}" $fcounts | awk '{print $14}')
	hostrna1=$(grep -P "repeat_region\t${hostnode}" $fcounts | awk '{print $15}')
	hostrna2=$(grep -P "repeat_region\t${hostnode}" $fcounts | awk '{print $16}')
	phagepro1=$(grep "${phagenode}" $fcounts | awk -F$'\t' '{print $13}' | paste -s -d+ | bc)
	phagepro2=$(grep "${phagenode}" $fcounts | awk -F$'\t' '{print $14}' | paste -s -d+ | bc)
	phagerna1=$(grep "${phagenode}" $fcounts | awk -F$'\t' '{print $15}' | paste -s -d+ | bc)
	phagerna2=$(grep "${phagenode}" $fcounts | awk -F$'\t' '{print $16}' | paste -s -d+ | bc)

	echo -e "${hostnode}\t${hostnodelen}\t${cstart}\t${cend}\t${phagenode}\t${phagenodelen}\t${pstart}\t${pend}\t${blastscore}\t${hostbin}\t${hostannot}\t${hostcov}\t${phagebin}\t${phageannot}\t${phagecov}\t${hostpro1}\t${hostpro2}\t${hostrna1}\t${hostrna2}\t${phagepro1}\t${phagepro2}\t${phagerna1}\t${phagerna2}" \
	>> ${NAME}_CRISPR_data.txt

done < $blastout
