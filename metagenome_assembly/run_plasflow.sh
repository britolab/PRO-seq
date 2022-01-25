#$ -S /bin/bash
#$ -N plasflow
#$ -V
#$ -o /workdir/users/acv46/mgmAssembly/log/plasflow_$JOB_ID.out
#$ -e /workdir/users/acv46/mgmAssembly/log/plasflow_$JOB_ID.err
#$ -wd /workdir/users/acv46/mgmAssembly
#$ -l h_vmem=100G
#$ -q long.q@cbsubrito
#$ -t 1

# last edited 13 Mar 2021, Albert Vill

WRK=/workdir/users/acv46/mgmAssembly/US3_10_CAB/metaSPAdes

##############
## PlasFlow ##
##############

echo "activating plasflow env"

export PATH=/programs/Anaconda2/bin:$PATH
export LD_LIBRARY_PATH=/programs/Anaconda2/lib:$LD_LIBRARY_PATH
source activate plasflow

#/programs/PlasFlow/scripts/filter_sequences_by_length.pl \
#	-input $WRK/contigs.fasta \
#	-output $WRK/contigs_1k.fasta \
#	-thresh 1000

echo "starting plasflow"
PlasFlow.py \
	--input $WRK/contigs_1k.fasta \
	--output $WRK/plasflow_0.99_out.tsv \
	--threshold 0.99
echo "plasflow finished"
#PlasFlow.py \
#        --input $WRK/contigs_1k.fasta \
#        --output $WRK/plasflow_0.7_out.tsv \
#        --threshold 0.7

source deactivate
