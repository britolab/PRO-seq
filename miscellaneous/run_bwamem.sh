#$ -S /bin/bash
#$ -N bwamem
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq2/log/bwamem_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq2/log/bwamem_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=60G
#$ -q long.q@cbsubrito
#$ -t 1-4

WRK=/workdir/users/acv46/stool_PROSeq2/deseq
BINDEX=/workdir/users/acv46/stool_PROSeq2/ref/haloferax/GCF_000337315.1_ASM33731v1_genomic.fna
FQ=/workdir/users/acv46/stool_PROSeq2/fastq/clean
export OMP_NUM_THREADS=6

DESIGN_FILE=$WRK/samples_PRO.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

cd $WRK

# get the sample name to find appropriate indices
samp=$(echo $NAME | awk -F_ '{print $1}')
# get new identifier to reduce name length
new=$(echo $NAME | sed 's/_rep//g' | awk -F_ '{print $1"_"$3}')

# do for HQ bin contig index

echo "skipping stuff"

: <<'END'

echo "shouldn't see this"

# Everything inside here is ignored
# easier than commenting out each line

bwa mem -a -v 2 \
	-t $OMP_NUM_THREADS \
	$BINDEX/${samp}_goodbins.short.fasta \
	${FQ}/${NAME}_R1.fastq.gz \
	${FQ}/${NAME}_R2.fastq.gz \
	> bam/${new}.goodbins_unsort.bam

samtools view -u bam/${new}.goodbins_unsort.bam | \
	samtools sort -@ $OMP_NUM_THREADS -o bam/${new}.goodbins.bam

rm bam/${new}.goodbins_unsort.bam

samtools index -@ $OMP_NUM_THREADS -b bam/${new}.goodbins.bam

# do for all contig index

bwa mem -a -v 2 \
        -t $OMP_NUM_THREADS \
        $BINDEX/${samp}_allcontigs.short.fasta \
        ${FQ}/${NAME}_R1.fastq.gz \
        ${FQ}/${NAME}_R2.fastq.gz \
        > bam/${new}.allcontigs_unsort.bam

samtools view -u bam/${new}.allcontigs_unsort.bam | \
        samtools sort -@ $OMP_NUM_THREADS -o bam/${new}.allcontigs.bam

rm bam/${new}.allcontigs_unsort.bam

samtools index -@ $OMP_NUM_THREADS -b bam/${new}.allcontigs.bam

END

# end ignored section

echo "aligning to Hfax reference"
# align PRO-seq reads to haloferax genome

bwa mem -a -v 2 \
        -t $OMP_NUM_THREADS \
        $BINDEX \
        ${FQ}/${NAME}_R1.fastq.gz \
        ${FQ}/${NAME}_R2.fastq.gz \
        > bam/${new}.hfax_unsort.bam

samtools view -u bam/${new}.hfax_unsort.bam | \
        samtools sort -@ $OMP_NUM_THREADS -o bam/${new}.hfax.bam

rm bam/${new}.hfax_unsort.bam

samtools index -@ $OMP_NUM_THREADS -b bam/${new}.hfax.bam
