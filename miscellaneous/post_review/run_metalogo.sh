#$ -S /bin/bash
#$ -N run_metalogo
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq3/log/metalogo_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq3/log/metalogo_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq3
#$ -l h_vmem=50G
#$ -q long.q@cbsubrito
#$ -t 1-2
#$ -tc 2
#$ -pe parenv 8

export OMP_NUM_THREADS=8
WRK=/workdir/users/acv46/stool_PROSeq3/figures/pause/metalogo
res=99
#fix=17

# peaks called with get_motifs R function
# peak results split to fastas with trunc_and_split R function
# see 3prime_pause.Rmd

DESIGN_FILE=${WRK}/jobs.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

fastas=${WRK}/${NAME}_fasta

mkdir -p $WRK/${NAME}_pdfout
mkdir -p $WRK/${NAME}_processed

source /home/acv46/miniconda3/bin/activate
conda activate metalogo

while read fasta; do

  bin=$(echo ${fasta} | sed "s/_3peaks.fasta//g" | sed "s/\./_/g")

  echo "running metalogo for ${NAME}, bin = ${bin}"

  metalogo \
    --seq_file ${fastas}/${fasta} \
    --seq_file_type fasta \
    --sequence_type dna \
    --task_name ${bin} \
    --group_resolution 0.${res} \
    --output_dir $WRK/${NAME}_pdfout \
    --fa_output_dir $WRK/${NAME}_processed \
    --output_name ${bin}_res${res} \
    --min_length 16 \
    --max_length 18 \
    --clustalo_bin /programs/clustalo \
    --fasttree_bin /home/acv46/fasttree/FastTree \
    --fasttreemp_bin /home/acv46/fasttree/FastTreeMP
    #--logo_format pdf

  echo "finished metalogo for ${NAME}, bin = ${bin}"

done < <(ls ${fastas} | grep "_3peaks.fasta")

