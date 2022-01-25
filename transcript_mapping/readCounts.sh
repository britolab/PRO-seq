#$ -S /bin/bash
#$ -N readcounts
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq2/log/readcounts_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq2/log/readcounts_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=10G
#$ -q long.q@cbsubrito
#$ -t 1-8
#$ -hold_jid 466363

#Set directories
WRK=/workdir/users/acv46/stool_PROSeq2
FQR=$WRK/fastq/raw
FQC=$WRK/fastq/clean

#Create design file of file names
DESIGN_FILE=$WRK/scripts/samples.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

raw=$FQR/${NAME}
bmtag=$FQC/${NAME}

# assumes raw files are gzipped, cleaned files are uncompressed
raw1=$(echo $(zcat ${raw}_R1.fastq|wc -l)/4|bc)
raw2=$(echo $(zcat ${raw}_R2.fastq|wc -l)/4|bc)
bmtag1=$(echo $(cat ${bmtag}_1.fastq|wc -l)/4|bc)
bmtag2=$(echo $(cat ${bmtag}_2.fastq|wc -l)/4|bc)

echo -e "${NAME}\t${raw1}\t${raw2}\t${bmtag1}\t${bmtag2}" \
        >> $FQC/readcounts_$JOB_ID.txt
