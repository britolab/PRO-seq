#$ -S /bin/bash
#$ -N combine_pileups
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq3/log/combine_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq3/log/combine_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq3
#$ -l h_vmem=100G
#$ -q long.q@cbsubrito
#$ -t 1
#$ -tc 1
#$ -pe parenv 8

# merge replicate bam files
# filter merged file
# create bedtools pileup

export OMP_NUM_THREADS=8

samp=US2
WRK=/workdir/users/acv46/stool_PROSeq3/transcriptomes/${samp}/coverage

echo "${samp} -- start"

cd $WRK

paste ${samp}_A_PROseq_fwd_depth.txt \
      ${samp}_A_PROseq_rev_depth.txt \
      ${samp}_A_TEXm_fwd_depth.txt \
      ${samp}_A_TEXm_rev_depth.txt \
      ${samp}_A_TEXp_fwd_depth.txt \
      ${samp}_A_TEXp_rev_depth.txt \
      ${samp}_B_PROseq_fwd_depth.txt \
      ${samp}_B_PROseq_rev_depth.txt \
      ${samp}_B_TEXm_fwd_depth.txt \
      ${samp}_B_TEXm_rev_depth.txt \
      ${samp}_B_TEXp_fwd_depth.txt \
      ${samp}_B_TEXp_rev_depth.txt \
      ${samp}_C_PROseq_fwd_depth.txt \
      ${samp}_C_PROseq_rev_depth.txt \
      ${samp}_C_TEXm_fwd_depth.txt \
      ${samp}_C_TEXm_rev_depth.txt \
      ${samp}_C_TEXp_fwd_depth.txt \
      ${samp}_C_TEXp_rev_depth.txt |
  awk '{print $1,$2,$3,$6,$9,$12,$15,$18,$21,$24,$27,$30,$33,$36,$39,$42,$45,$48,$51,$54}' |
  tr [[:blank:]] "\t" | \
  sed "s/_length[^\t]*//g" | sed "s/NODE_//g"  | \
  sed '1i contig\tposition\tPRO_A_f\tPRO_A_r\tTEXm_A_f\tTEXm_A_r\tTEXp_A_f\tTEXp_A_r\tPRO_B_f\tPRO_B_r\tTEXm_B_f\tTEXm_B_r\tTEXp_B_f\tTEXp_B_r\tPRO_C_f\tPRO_C_r\tTEXm_C_f\tTEXm_C_r\tTEXp_C_f\tTEXp_C_r' \
  > ${samp}_combined_depth.txt

echo "${samp} -- done!"

