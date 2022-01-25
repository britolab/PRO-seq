#$ -S /bin/bash
#$ -N anvio
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/anvio_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/anvio_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=80G
#$ -q long.q@cbsubrito
#$ -t 1-5
#$ -hold_jid 468001

# this script takes in metaSPAdes assemblies (contigs.fasta)
# and reads mapped to those assemblies (sample.bam)
# and creates anvio databases for visualization

# https://merenlab.org/2016/06/22/anvio-tutorial-v2/
# https://merenlab.org/2015/06/10/combining-omics-data/

# once contig and profile databases are made, visualize with the following commands
## ssh -L 8080:localhost:8080 acv46@cbsubrito.tc.cornell.edu
## conda activate /workdir/acv46/anvio7
## anvi-interactive -p SAMPLE/PROFILE.db -c contigs.db
# then open connection in chrome at localhost:8080

WRK=/workdir/users/acv46/stool_PROSeq2/anvio/US2_prevo

DESIGN_FILE=$WRK/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

new=$(echo $NAME | sed "s/_rep//g")

ANV=$WRK/${new}
CTG=/workdir/users/acv46/stool_PROSeq2/ref/US2_bins/merged/prevo_bins.fa

# start anvio
source $HOME/miniconda3/bin/activate
conda activate /workdir/acv46/anvio7
export OMP_NUM_THREADS=4

# check for illegal characters in fasta deflines

cd $WRK

# generate contigs database
#anvi-gen-contigs-database \
#	-f $CTG \
#	-o contigs.db \
#	-n 'US2 Prevotella contigs database' \
#	-T $OMP_NUM_THREADS

# run hmms
#anvi-run-hmms \
#	-c contigs.db \
#	-T $OMP_NUM_THREADS

# annotate with COGS
# run anvi-setup-ncbi-cogs if first time
#anvi-run-ncbi-cogs \
#	-c contigs.db \
#       -T $OMP_NUM_THREADS

# sort and index bam files
# create anvio profiles

anvi-profile \
	-i bam/${new}.bam \
	-c contigs.db \
	-M 2000 \
	-T $OMP_NUM_THREADS \
	-o ${new} \
	-S ${new}

# merge without clustering for transcriptomics data

#anvi-merge ./profile_*/PROFILE.db -c contigs.db -o MERGED \
#           --skip-hierarchical-clustering \
#           --skip-concoct-binning

conda deactivate
