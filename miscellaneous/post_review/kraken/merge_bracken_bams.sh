#$ -S /bin/bash
#$ -N merge_bracken
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq3/log/merge_bracken_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq3/log/merge_bracken_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq3
#$ -l h_vmem=50G
#$ -q long.q@cbsubrito
#$ -t 1
#$ -tc 1
#$ -M acv46@cornell.edu
#$ -pe parenv 8

# this script runs combine_bracken_outputs.py
# using the output of bracken.sh

WRK=/workdir/users/acv46/stool_PROSeq3/kraken
export OMP_NUM_THREADS=8
# see /workdir/refdbs/kraken/DB/_README for database install details

source /home/acv46/miniconda3/bin/activate
conda activate kraken2

cd ${WRK}

levels=P,C,O,F,G,S,S1

for level in $(echo $levels | sed "s/,/ /g"); do

	combine_bracken_outputs.py \
	--files ./kraken2_bam/*/*.bracken_${level}.txt \
	--output ./merged_bracken_bams/bracken_merged_${level}.txt \
	--names US2_A_PROseq,US2_A_RNAseq,US2_A_TEXm,US2_A_TEXp,US2_B_PROseq,US2_B_RNAseq,US2_B_TEXm,US2_B_TEXp,US2_C_PROseq,US2_C_TEXm,US2_C_TEXp,US2_mgm,US3_A_PROseq,US3_A_RNAseq,US3_A_TEXm,US3_A_TEXp,US3_B_PROseq,US3_B_RNAseq,US3_B_TEXm,US3_B_TEXp,US3_C_PROseq,US3_C_TEXm,US3_C_TEXp,US3_mgm

done

conda deactivate
