#$ -S /bin/bash
#$ -N vs2_train
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq2/log/vs2train_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq2/log/vs2train_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=100G
#$ -q long.q@cbsubrito
#$ -t 1

# installation
# conda install virsorter
# conda update virsorter
# virsorter setup -d db -j 4

# this script runs the virsorter train-feature snakemake pipeline
# goal: create classifier for human Gut Phage Database

WRK=/workdir/users/acv46/stool_PROSeq2/phage/virsorter
VSORT=/home/acv46/miniconda3/bin/virsorter
#export PATH=/programs/prodigal-2.6.3:$PATH
export OMP_NUM_THREADS=16

#DESIGN_FILE=$WRK/assemblies.txt
#        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
#        NAME=`basename ${DESIGN}`

fasta=/workdir/refdbs/GutPhageDatabase/GPD_sequences.fa

#echo "generating prodigal training file for RBS identification"

#prodigal \
#        -i $fasta \
#	-g 11 \
#        -t $WRK/GPD_prodigal.train

source /home/acv46/miniconda3/bin/activate

cd $WRK

echo "generating virsorter feature file for Gut Phage Database"

# default hmms
# no hallmark genes
# no prodigal RBS motif training model

$VSORT train-feature \
	--seqfile /workdir/refdbs/GutPhageDatabase/GPD_sequences.fa \
	--hmm db/hmm/viral/combined.hmm \
	--frags-per-genome 5 \
	--jobs 16 \
	-w GPD-feature.out \
	--max-orf-per-seq 40



