#$ -S /bin/bash
#$ -N mapHiC
#$ -V
#$ -o /workdir/users/acv46/mgmAssembly/log/mapHiC_$JOB_ID.out
#$ -e /workdir/users/acv46/mgmAssembly/log/mapHiC_$JOB_ID.err
#$ -wd /workdir/users/acv46/mgmAssembly
#$ -l h_vmem=80G
#$ -q long.q@cbsubrito
#$ -t 1-5

# last edited 29 Mar 2021, Albert Vill

WRK=/workdir/users/acv46/mgmAssembly

#################################
## Switches and Global Options #############################################
#################################                                          #
									   #
# increase thread restriction of cluster				   #
## default = 1								   #
## this value is used by metaSPAdes, BWA, CONCOCT, metaBAT,		   #
### MaxBin, DASTool, checkM, and GTDB-Tk 				   #
export OMP_NUM_THREADS=8						   #

# to receive emails when the script fails or finishes,			   #
# add email adress here							   #
# otherwise, leave as email=NULL					   #
## NOTE -- sends one email per task					   #
email=acv46@cornell.edu							   #
									   #
# location of jobtracker.sh script					   #
jobtracker=/workdir/users/acv46/mgmAssembly/scripts/jobtracker.sh	   #
									   #
############################################################################

DESIGN_FILE=$WRK/names.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

SAMP=$WRK/${NAME}_CAB
mkdir -p $SAMP

FASTQ1_mgm=$SAMP/fastq/${NAME}.clean_1.fastq
FASTQ2_mgm=$SAMP/fastq/${NAME}.clean_2.fastq
FASTQ1_D=$SAMP/D_proxlig/${NAME}_D.clean_1.fastq
FASTQ2_D=$SAMP/D_proxlig/${NAME}_D.clean_2.fastq
FASTQ1_B=$SAMP/B_proxlig/${NAME}_B.clean_1.fastq
FASTQ2_B=$SAMP/B_proxlig/${NAME}_B.clean_2.fastq
FASTQ1_N=$SAMP/N_proxlig/${NAME}_N.clean_1.fastq
FASTQ2_N=$SAMP/N_proxlig/${NAME}_N.clean_2.fastq

echo "${NAME} -- running plasflow on metagenomic contigs"

MSPA=$SAMP/metaSPAdes
mkdir -p $MSPA
cd $MSPA

cont1=$(grep ">" $MSPA/contigs.fasta | wc -l)

echo "${NAME} -- starting with ${cont1} contigs"

export PATH=/programs/Anaconda2/bin:$PATH
export LD_LIBRARY_PATH=/programs/Anaconda2/lib:$LD_LIBRARY_PATH
source activate plasflow

/programs/PlasFlow/scripts/filter_sequences_by_length.pl \
	-input $MSPA/contigs.fasta \
	-output $MSPA/plasflow_contigs.fasta \
	-thresh 1000

cont2=$(grep ">" $MSPA/plasflow_contigs.fasta | wc -l)

echo "${NAME} -- ${cont2} contigs > 1kb"

thresh=0.95

PlasFlow.py \
	--input $MSPA/plasflow_contigs.fasta \
	--output $MSPA/plasflow_0.95_out.tsv \
	--threshold $thresh

source deactivate

cont_chr=$(grep ">" $MSPA/plasflow_${thresh}_out.tsv_chromosomes.fasta | wc -l)
cont_pls=$(grep ">" $MSPA/plasflow_${thresh}_out.tsv_plasmids.fasta | wc -l)
cont_unc=$(grep ">" $MSPA/plasflow_${thresh}_out.tsv_unclassified.fasta | wc -l)

sed -i "s/>/>CHROM_/g" $MSPA/plasflow_${thresh}_out.tsv_chromosomes.fasta
sed -i "s/>/>PLAS_/g" $MSPA/plasflow_${thresh}_out.tsv_plasmids.fasta
sed -i "s/>/>UNCLASS_/g" $MSPA/plasflow_${thresh}_out.tsv_unclassified.fasta

cat $MSPA/plasflow_${thresh}_out.tsv_chromosomes.fasta \
	$MSPA/plasflow_${thresh}_out.tsv_plasmids.fasta \
	$MSPA/plasflow_${thresh}_out.tsv_unclassified.fasta \
	> $MSPA/plasflow_${thresh}_labeled.fasta

echo -e "${NAME} -- plasflow results for threshold=${thresh}\n----> ${cont_chr} chromosomal contigs\n----> ${cont_pls} plasmid contigs\n----> ${cont_unc} unclassified contigs"

#####

MAP=$SAMP/map
mkdir -p $MAP
cd $MAP

export PATH=/programs/bwa-0.7.17:$PATH
export PATH=/programs/samtools-1.11/bin:$PATH
DEPTH=/programs/metabat/jgi_summarize_bam_contig_depths

echo "${NAME} -- making bwa index of contigs"

bwa index -a bwtsw -p plasflow_contigs $MSPA/plasflow_${thresh}_labeled.fasta

#####

echo "${NAME} -- mapping metagenome reads"

bwa mem -a -v 2 -t $OMP_NUM_THREADS $MAP/plasflow_contigs \
        $FASTQ1_mgm $FASTQ2_mgm \
        > $MAP/mgm.unsort.bam
samtools view -u $MAP/mgm.unsort.bam | \
        samtools sort -@ 4 -o $MAP/mgm.bam

rm $MAP/mgm.unsort.bam

samtools index -@ 4 -b $MAP/mgm.bam

DEPTH=/programs/metabat/jgi_summarize_bam_contig_depths
$DEPTH --outputDepth $MAP/mgm.depth.txt $MAP/mgm.bam

echo "${NAME} -- mapping finished, metagenome"

#####

echo "${NAME} -- mapping D reads"

bwa mem -a -v 2 -t $OMP_NUM_THREADS $MAP/plasflow_contigs \
        $FASTQ1_D $FASTQ2_D \
        > $MAP/D.unsort.bam
samtools view -u $MAP/D.unsort.bam | \
        samtools sort -@ 4 -o $MAP/D.bam

rm $MAP/D.unsort.bam

samtools index -@ 4 -b $MAP/D.bam

DEPTH=/programs/metabat/jgi_summarize_bam_contig_depths
$DEPTH --outputDepth $MAP/D.depth.txt $MAP/D.bam

echo "${NAME} -- mapping finished, D"

#####

echo "${NAME} -- mapping B reads"

bwa mem -a -v 2 -t $OMP_NUM_THREADS $MAP/plasflow_contigs \
        $FASTQ1_B $FASTQ2_B \
        > $MAP/B.unsort.bam
samtools view -u $MAP/B.unsort.bam | \
        samtools sort -@ 4 -o $MAP/B.bam

rm $MAP/B.unsort.bam

samtools index -@ 4 -b $MAP/B.bam

DEPTH=/programs/metabat/jgi_summarize_bam_contig_depths
$DEPTH --outputDepth $MAP/B.depth.txt $MAP/B.bam

echo "${NAME} -- mapping finished, B"

#####

echo "${NAME} -- mapping N reads"

bwa mem -a -v 2 -t $OMP_NUM_THREADS $MAP/plasflow_contigs \
        $FASTQ1_N $FASTQ2_N \
        > $MAP/N.unsort.bam
samtools view -u $MAP/N.unsort.bam | \
        samtools sort -@ 4 -o $MAP/N.bam

rm $MAP/N.unsort.bam

samtools index -@ 4 -b $MAP/N.bam

DEPTH=/programs/metabat/jgi_summarize_bam_contig_depths
$DEPTH --outputDepth $MAP/N.depth.txt $MAP/N.bam

echo "${NAME} -- mapping finished, N"

#####

echo "${NAME} -- starting metaphlan"

MPHL=$SAMP/metaphlan
mkdir -p $MPHL
cd $MPHL

export PATH=/programs/MetaPhlAn-2.0:/programs/MetaPhlAn-2.0/utils:$PATH

metaphlan2.py \
	--input_type fastq \
	--bowtie2db /workdir/users/fnn3/references/db_v20/mpa_v20_m200 \
	--bowtie2out ${NAME}_mgm.bowtie2.bz2 \
	--nproc $OMP_NUM_THREADS \
	$FASTQ1_mgm,$FASTQ2_mgm \
	${NAME}_mgm_profile.txt

metaphlan2.py \
        --input_type fastq \
        --bowtie2db /workdir/users/fnn3/references/db_v20/mpa_v20_m200 \
        --bowtie2out ${NAME}_D.bowtie2.bz2 \
        --nproc $OMP_NUM_THREADS \
        $FASTQ1_D,$FASTQ2_D \
        ${NAME}_D_profile.txt

metaphlan2.py \
        --input_type fastq \
        --bowtie2db /workdir/users/fnn3/references/db_v20/mpa_v20_m200 \
        --bowtie2out ${NAME}_B.bowtie2.bz2 \
        --nproc $OMP_NUM_THREADS \
        $FASTQ1_B,$FASTQ2_B \
        ${NAME}_B_profile.txt

metaphlan2.py \
        --input_type fastq \
        --bowtie2db /workdir/users/fnn3/references/db_v20/mpa_v20_m200 \
        --bowtie2out ${NAME}_N.bowtie2.bz2 \
        --nproc $OMP_NUM_THREADS \
        $FASTQ1_N,$FASTQ2_N \
        ${NAME}_N_profile.txt
