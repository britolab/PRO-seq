#$ -S /bin/bash
#$ -N mnmE
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq2/log/mnmE_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq2/log/mnmE_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=10G
#$ -q long.q@cbsubrito
#$ -t 1-11

# see /workdir/users/acv46/stool_PROSeq2/deseq/featureCounts/US2_goodbins/_README

WRK=/workdir/users/acv46/stool_PROSeq2/deseq/featureCounts/US2_goodbins
export OMP_NUM_THREADS=6

DESIGN_FILE=$WRK/goodbins_90_10.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

cd $WRK

# get contigs for bin

grep -P "^${NAME}\t" US2_bin2contig.txt | \
	awk '{print $2}' > bin_${NAME}_contiglist.tmp

# get counts across contigs for each bam
while read contig; do
	p1=$(grep "${contig}" US2_PRO1.primarymapped.txt | awk '{print $1}')
	p2=$(grep "${contig}" US2_PRO2.primarymapped.txt | awk '{print $1}')
	r1=$(grep "${contig}" US2_RNA1.primarymapped.txt | awk '{print $1}')
        r2=$(grep "${contig}" US2_RNA2.primarymapped.txt | awk '{print $1}')
	echo -e "${NAME}.fa\t${contig}\t${p1}\t${p2}\t${r1}\t${r2}" >> bin_contig_counts.txt
done < bin_${NAME}_contiglist.tmp

rm bin_${NAME}_contiglist.tmp


