#$ -S /bin/bash
#$ -N anvio_profile
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/anvio_profile_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/anvio_profile_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=50G
#$ -q long.q@cbsubrito
#$ -t 1-5

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

WRK=/workdir/users/acv46/stool_PROSeq2/anvio/US3_3Nov2020

DESIGN_FILE=$WRK/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

source $HOME/miniconda3/bin/activate
conda activate /workdir/acv46/anvio7
export OMP_NUM_THREADS=4

cd $WRK

samp=$(echo $NAME | sed "s/US3_3Nov2020_//g" | sed "s/_rep//g")

anvi-profile \
	-i ${NAME}.sort.bam \
	-c contigs.db \
	-M 2000 \
	-T $OMP_NUM_THREADS \
	-o $WRK/${samp} \
	-S ${samp}

conda deactivate
