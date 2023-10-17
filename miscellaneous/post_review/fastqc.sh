#$ -S /bin/bash
#$ -N fastqc
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq3/log/fastqc_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq3/log/fastqc_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq3/
#$ -l h_vmem=20G
#$ -q long.q@cbsubrito
#$ -t 1-22
#$ -tc 10
#$ -pe parenv 4

# last edited 3 Nov 2022, Albert Vill

# this scripts runs fastqc on transcriptomics reads
# RNA-seq, dRNA-seq (TEX-treated), PRO-seq

WRK=/workdir/users/acv46/stool_PROSeq3/transcriptomes
FQ=/home/britolab/data/NextSeq_10458816
export PATH=/programs/FastQC-0.11.8:$PATH
export OMP_NUM_THREADS=4

# Create design file of file names
DESIGN_FILE=$WRK/names.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

OUT=$WRK/fastqc
mkdir -p $OUT
cd $OUT

raw1=$(grep -P "\t${NAME}" $WRK/names_to_reads.txt | cut -f1 | head -n1 | sed "s/.fastq.gz//g")
raw2=$(grep -P "\t${NAME}" $WRK/names_to_reads.txt | cut -f1 | tail -n1 | sed "s/.fastq.gz//g")

fastqc \
	-o $OUT \
	--noextract \
	-t ${OMP_NUM_THREADS} \
	-a $WRK/adapters.txt \
	${FQ}/${raw1}.fastq.gz \
	${FQ}/${raw2}.fastq.gz

mv ${raw1}_fastqc.html ${NAME}_R1_fastqc.html
mv ${raw2}_fastqc.html ${NAME}_R2_fastqc.html

rm ${raw1}_fastqc.zip
rm ${raw2}_fastqc.zip

