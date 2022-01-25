#$ -S /bin/bash
#$ -N roary
#$ -V
#$ -o /workdir/users/acv46/ben_May2021/log/roary_$JOB_ID.out
#$ -e /workdir/users/acv46/ben_May2021/log/roary_$JOB_ID.err
#$ -wd /workdir/users/acv46/ben_May2021
#$ -l h_vmem=60G
#$ -q long.q@cbsubrito2
#$ -t 1-292

# prokka installed via conda

# roary installed via conda
## source /home/acv46/miniconda3/bin/activate
## conda config --add channels r
## conda config --add channels defaults
## conda config --add channels conda-forge
## conda config --add channels bioconda
## conda create --name roary
## conda activate roary
## conda install roary

WRK=/workdir/users/acv46/ben_May2021/vamb/small_NN
BIN=$WRK/bins
KRA=$WRK/krakenuniq
export OMP_NUM_THREADS=8

# cd bins
# ls ./*.fna | \
# awk -F$"C" '{print $2}' | \
# sort | uniq \
# > ../mergedbins.txt

DESIGN_FILE=$WRK/mergedbins.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

genome=$(echo "C${NAME}" | sed "s/.fna//g")
OUT=$WRK/roary/${genome}
mkdir -p $OUT
ls $BIN | grep "S.*${genome}.fna$" > $OUT/bins.txt

echo "${NAME} -- this identifier represents $(cat $OUT/bins.txt | wc -l) bin(s)"

# get kraken genus of bin set
# run prokka using genus-specific database if all bins are same genus

genus=$(grep -m 1 "genus" $KRA/S*${genome}.fna.report.txt | awk '{print $NF}' | uniq)
numgen=$(grep -m 1 "genus" $KRA/S*${genome}.fna.report.txt | awk '{print $NF}' | uniq | wc -l)

source /home/acv46/miniconda3/bin/activate

if [[ "$numgen" -eq 1 ]]; then

	echo "${NAME} -- all bins annotated as ${genus}"

	while read fasta; do

		fname=$(echo $fasta | sed "s/.fna//g")

		sed "s/_length.*$//g" $BIN/$fasta \
		> $OUT/${fname}_edit.fna
		prokka \
			--kingdom bacteria \
			--outdir $OUT/temp_${fname} \
			--genus $genus \
			--locustag $fname \
			--cpus $OMP_NUM_THREADS \
			$OUT/${fname}_edit.fna

		mkdir -p $OUT/gff
		mv $OUT/temp_${fname}/*.gff $OUT/gff/${fname}.gff

		rm $OUT/${fname}_edit.fna
		rm -r $OUT/temp_${fname}

	done < $OUT/bins.txt

else

	echo "${NAME} -- bins have different kraken annotaions, running prokka without --genus flag"

 	while read fasta; do

                fname=$(echo $fasta | sed "s/.fna//g")

                sed "s/_length.*$//g" $BIN/$fasta \
                > $OUT/${fname}_edit.fna
                prokka \
                        --kingdom bacteria \
                        --outdir $OUT/temp_${fname} \
                        --locustag $fname \
                        --cpus $OMP_NUM_THREADS \
                        $OUT/${fname}_edit.fna

                mkdir -p $OUT/gff
                mv $OUT/temp_${fname}/*.gff $OUT/gff/${fname}.gff

                rm $OUT/${fname}_edit.fna
		rm -r $OUT/temp_${fname}

        done < $OUT/bins.txt

fi

export PATH=/programs/R-4.0.0/bin:$PATH
export PATH=/home/acv46/R/x86_64-pc-linux-gnu-library/4.0:$PATH

echo "${NAME} -- starting roary"

numgff=$(ls $OUT/gff/*.gff | wc -l)

if [[ "$numgff" -eq 1 ]]; then

	echo "${NAME} -- only one gff, need >=2 to run roary"
	exit 0

else

	echo "${NAME} -- starting roary using ${numgff} bin annotations"

	conda activate roary

	# enter a tmpdir to avoid file locking
	mkdir -p $WRK/tmpdir
	cd $WRK/tmpdir

	roary \
		-p $OMP_NUM_THREADS \
		-f $OUT/out \
		-e \
		-n \
		-i 95 \
		-r \
		$OUT/gff/*.gff

	echo "${NAME} -- roary finished"

	conda deactivate

fi
