#$ -S /bin/bash
#$ -N krakenuniq
#$ -V
#$ -o /workdir/users/acv46/ben_May2021/log/krakenuniq_$JOB_ID.out
#$ -e /workdir/users/acv46/ben_May2021/log/krakenuniq_$JOB_ID.err
#$ -wd /workdir/users/acv46/ben_May2021
#$ -l h_vmem=300G
#$ -q long.q@cbsubrito2
#$ -t 1-3

# holy shit kraken uses a lot of memory
# allocate at least 290G for DB_subset

# krakenuniq installed via conda
## source /home/acv46/miniconda3/bin/activate
## conda create --name krakenuniq
## conda activate krakenuniq
## conda install krakenuniq

WRK=/workdir/users/acv46/ben_May2021/vamb
export OMP_NUM_THREADS=32
kdb=/workdir/users/agk85/tools/krakenuniq/DB_subset

source /home/acv46/miniconda3/bin/activate
conda activate krakenuniq

echo "loading database into memory"

krakenuniq \
	--db $kdb \
	--preload \
	--threads $OMP_NUM_THREADS

DESIGN_FILE=$WRK/models.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

bins=$WRK/${NAME}/bins

OUT=$WRK/${NAME}/krakenuniq
mkdir -p $OUT
cd $OUT

echo "${NAME} -- starting krakenuniq"

ls $bins > $OUT/bins.txt
bincount=$(cat $OUT/bins.txt | wc -l)

while read bin; do

	echo "${NAME}: running krakenuniq for bin $(grep -wn "${bin}" $OUT/bins.txt | awk -F ":" '{print $1}') of ${bincount} -- ${bin}"

	krakenuniq \
		--report-file ${bin}.report.txt \
		--db $kdb \
		--threads $OMP_NUM_THREADS \
		--fasta-input \
		--output ${bin}.out.txt \
		${bins}/${bin}

done < $OUT/bins.txt

echo "${NAME} -- finished krakenuniq"

conda deactivate
