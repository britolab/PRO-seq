#$ -S /bin/bash
#$ -N bakta
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq3/log/bakta_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq3/log/bakta_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq3
#$ -l h_vmem=100G
#$ -q long.q@cbsubrito
#$ -t 1-149
#$ -tc 4
#$ -pe parenv 8

# annotate MAGs with bakta
# give at least 50G to run diamond

WRK=/workdir/users/acv46/stool_PROSeq3
export OMP_NUM_THREADS=8

# generate file of absolute bin paths
## cd /workdir/users/acv46/stool_PROSeq3/assembly
## find ${PWD} -name "*.fa" | grep "DASTool_bins" > binlist.txt

LIST=$WRK/assembly/binlist.txt
  DESIGN=$(sed -n "${SGE_TASK_ID}p" $LIST)
  NAME="$DESIGN"

# debugging exit
## echo ${SGE_TASK_ID}
## echo ${NAME}
## exit 1

samp=$(basename $(dirname ${NAME}) | cut -f1 -d'_')
bin=$(basename ${NAME} | sed "s/.fa//g")
stats=$WRK/assembly/${samp}_25May2022/DASTool/${samp}_25May2022_BINSTATS.txt
OUT=$WRK/annotation/${samp}_bins_bakta_1.7/${bin}
mkdir -p ${OUT}
mkdir -p ${OUT}/temp
cd ${OUT}

# get taxa for bin
genus=$(grep -P "^${bin}\t" ${stats} | awk -F$'\t' '{print $15}')
species=$(grep -P "^${bin}\t" ${stats} | awk -F$'\t' '{print $16}' | sed 's/^.* //')

# bashrc alias
source /home/acv46/miniconda3/bin/activate
conda activate bakta

echo "START: Task ${SGE_TASK_ID}/149, running bakta on sample ${samp} bin ${bin}"

bakta \
  -d /workdir/refdbs/bakta/db \
  -p ${bin} \
  -o ${OUT} \
  -m 1000 \
  --genus ${genus} \
  --species ${species} \
  -t ${OMP_NUM_THREADS} \
  --keep-contig-headers \
  --tmp-dir ${OUT}/temp \
  ${NAME}

rm -rf ${OUT}/temp
echo "FINISHED: Task ${SGE_TASK_ID}/149, running bakta on sample ${samp} bin ${bin}"
