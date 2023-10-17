#$ -S /bin/bash
#$ -N bracken
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq3/log/bracken_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq3/log/bracken_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq3
#$ -l h_vmem=20G
#$ -q long.q@cbsubrito
#$ -t 1-10
#$ -tc 5
#$ -M acv46@cornell.edu
#$ -pe parenv 8

# runs bracken on kraken2 output to get relative abundance measures
# see kraken2.sh

# kraken2 installed via conda
## source /home/acv46/miniconda3/bin/activate
## conda create --name kraken2
## conda activate kraken2
## conda install -c bioconda kraken2

WRK=/workdir/users/acv46/stool_PROSeq3/kraken
export OMP_NUM_THREADS=8
# see /workdir/refdbs/kraken/HumGut_DB/_README for database install details
kdb=/workdir/refdbs/kraken2db
kout=$WRK/kraken2_old_bam

DESIGN_FILE=$WRK/old_names.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

OUT=$kout/${NAME}
mkdir -p $OUT
cd $OUT

echo "${NAME} -- starting bracken"

readlength=$(awk 'NR%4==2 {print length}' ${NAME}_converted.fastq | sort -n | uniq -c | sort -rh | head -1 | cut -f2)

source /home/acv46/miniconda3/bin/activate
conda activate kraken2

levels=P,C,O,F,G,S,S1

for level in $(echo $levels | sed "s/,/ /g"); do

	echo "--> ${NAME} -- running bracken for taxonomic level ${level}"

	bracken \
		-d $kdb \
        	-i ${NAME}.report.txt \
        	-o ${NAME}.bracken_${level}.txt \
        	-r 100 \
        	-l ${level}

done

conda deactivate

echo "${NAME} -- finished bracken, script done"
