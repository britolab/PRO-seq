#$ -S /bin/bash
#$ -N combine_pileups
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq3/log/combine_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq3/log/combine_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq3
#$ -l h_vmem=150G
#$ -q long.q@cbsubrito
#$ -t 1
#$ -tc 1
#$ -pe parenv 8

# merge replicate bam files
# filter merged file
# create bedtools pileup

export OMP_NUM_THREADS=8
WRK=/workdir/users/acv46/stool_PROSeq3/transcriptomes/pileup

qscore=30

List=$WRK/../base_US2.txt
  DESIGN=$(sed -n "${SGE_TASK_ID}p" $List)
  NAME=`basename "$DESIGN"`

cd $WRK

echo "${NAME} -- start combine"

paste ${NAME}_PROseq_plus_q${qscore}_full.cov ${NAME}_PROseq_minus_q${qscore}_full.cov \
      ${NAME}_PROseq_plus_q${qscore}_3p.cov ${NAME}_PROseq_minus_q${qscore}_3p.cov \
      ${NAME}_PROseq_plus_q${qscore}_5p.cov ${NAME}_PROseq_minus_q${qscore}_5p.cov \
      ${NAME}_RNAseq_plus_q${qscore}_full.cov ${NAME}_RNAseq_minus_q${qscore}_full.cov \
      ${NAME}_TEXp_plus_q${qscore}_full.cov ${NAME}_TEXp_minus_q${qscore}_full.cov \
      ${NAME}_TEXm_plus_q${qscore}_full.cov ${NAME}_TEXm_minus_q${qscore}_full.cov |
  awk '{print $1,$2,$3,$6,$9,$12,$15,$18,$21,$24,$27,$30,$33,$36}' | \
  tr [[:blank:]] "\t" | \
  sed "s/_length[^\t]*//g" | sed "s/NODE_//g"  | \
  sed '1i contig\tposition\tPRO_plus_full\tPRO_minus_full\tPRO_plus_3\tPRO_minus_3\tPRO_plus_5\tPRO_minus_5\tRNA_plus\tRNA_minus\tTEXp_plus\tTEXp_minus\tTEXm_plus\tTEXm_minus' \
  > ${NAME}_combined_q${qscore}_correct.txt

echo "${NAME} -- done!"

# remove intermediate files
rm ${NAME}_*_q${qscore}_*.cov
rm ${NAME}_*_q${qscore}.bam
