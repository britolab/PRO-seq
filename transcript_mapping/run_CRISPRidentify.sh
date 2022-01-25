#$ -S /bin/bash
#$ -N CRISPRid
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq2/log/CRISPRid_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq2/log/CRISPRid_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=80G
#$ -q long.q@cbsubrito
#$ -t 1-2

## Doesn't work, too many bugs, too long to run

# runs CRISPRidentify tool
# installed via conda
## cd /home/britolab/acv46/CRISPRidentify
## wget https://raw.githubusercontent.com/BackofenLab/CRISPRidentify/master/environment.yml
## source /home/acv46/miniconda3/bin/activate
## conda env create -f environment.yml
## conda activate crispr_identify_env

WRK=/workdir/users/acv46/stool_PROSeq2
#BLASTDB=$WRK/phage/virsorter
CID=/home/britolab/acv46/CRISPRidentify/CRISPRidentify-master
OUT=$WRK/phage/CRISPRidentify_out
CONTIG=$WRK/deseq/ref

export OMP_NUM_THREADS=20

DESIGN_FILE=$OUT/../allcontigs.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

source /home/acv46/miniconda3/bin/activate
conda activate crispr_identify_env

cd $CID

#OUT1=$CID/out/${NAME}
#mkdir -p $OUT1

python $CID/CRISPRidentify.py \
	--file $CID/${NAME} \
	--cpu $OMP_NUM_THREADS \
	--cas True \
	--is_element True \
	#--result_folder $OUT \
	--fasta_report True

conda deactivate

#mv $CID/${NAME} $OUT/


