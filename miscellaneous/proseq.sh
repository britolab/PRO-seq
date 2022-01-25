#$ -S /bin/bash
#$ -N proseq2.0
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/proseq_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/proseq_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2/fastq/clean
#$ -l h_vmem=180G
#$ -q long.q@cbsubrito
#$ -t 1-4

#The workdir above (-wd) specifies the location of the input reads
#Input paired-end reads must end in _R1.fastq.gz, _R2.fastq.gz

WRK=/workdir/users/acv46/stool_PROSeq2
JOB=$WRK/scripts/samples.txt
PRO=/workdir/users/acv46/stool_PROSeq/proseq2.0/proseq2.0-master
OUT=$WRK/out3
TMP=$OUT/temp
FASTQ=$WRK/fastq/clean

mkdir -p $OUT
mkdir -p $TMP

#generate bwaIndex using "bwa index" (generated automatically during metaSPADES assembly)
#generate chromInfo file using Get_Lengths.sh

export PATH=/programs/seqtk:$PATH
export PATH=/programs/bin/fastx:$PATH
export PATH=/programs/prinseq-lite-0.20.2:$PATH
export PATH=/programs/bwa-0.7.17:$PATH
export PATH=/programs/samtools-1.9/bin:$PATH
export PATH=/programs/bedtools2-2.28.0/bin:$PATH
export PATH=/programs/bedops-2.4.35/bin:$PATH
export PATH=/programs/kentUtils/bin:$PATH
export OMP_NUM_THREADS=8

LIST=$JOB
  DESIGN=$(sed -n "${SGE_TASK_ID}p" $LIST)
  NAME=`basename "$DESIGN"`

echo "${NAME} -- starting proseq2.0 pipeline"

mkdir -p $OUT/${NAME}

# access bwa index

echo "${NAME} -- locating existing BWA index"

base=$(echo ${NAME} | cut -f1,2 -d'_')
MGM=$OUT/${base}_ref
bwaIndex=$MGM/${base}_bins_plasmids

# create chromInfo file

echo "${NAME} -- creating chromInfo file for metaSPAdes assembly"

contigs=$MGM/${base}_bins_plasmids.fa
getlengths=$WRK/scripts/Get_Lengths.sh
chromInfo=${bwaIndex}.chromInfo
sh $getlengths $contigs $chromInfo

# process and map reads

mv $FASTQ/${NAME}_1.fastq $FASTQ/${NAME}_R1.fastq
mv $FASTQ/${NAME}_2.fastq $FASTQ/${NAME}_R2.fastq

READ1=$FASTQ/${NAME}_R1.fastq
READ2=$FASTQ/${NAME}_R2.fastq

gzip $READ1
gzip $READ2

if [[ "$NAME" == *"PRO"* ]]; then

	echo "${NAME} -- starting proseq2.0 for PRO-seq reads"

	cd $FASTQ

	bash $PRO/proseq2.0_edited.bsh \
		-i $bwaIndex \
		-c $chromInfo \
		-PE \
		--RNA3=R1_5prime \
		--ADAPT1=GATCGTCGGACTGTAGAACTCTGAAC \
		--ADAPT2=CCTTGGCACCCGAGAATTCCA \
		--UMI1=6 \
		--ADD_B1=0 \
		--opposite-strand=TRUE \
		-T $TMP \
		-O $OUT/${NAME} \
		-I ${NAME} \
		--thread=${OMP_NUM_TREADS}

elif [[ "$NAME" == *"RNA"* ]]; then

	echo "${NAME} -- starting proseq2.0 for RNA-seq reads"

	cd $FASTQ

        bash $PRO/proseq2.0_edited.bsh \
                -i $bwaIndex \
                -c $chromInfo \
                -PE \
                --RNA5=R1_5prime \
		--ADAPT2=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
		--ADAPT1=GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT \
                --Force_deduplicate=FALSE \
                --opposite-strand=FALSE \
                -T $TMP \
                -O $OUT/${NAME} \
                -I ${NAME} \
                --thread=${OMP_NUM_THREADS}

fi

echo "${NAME} -- proseq2.0 pipeline done, check that mapping is correct orientation"

#for PE PRO-seq, set --RNA3=R1_5prime
#for PE PRO-seq workflow, include
	#--UMI1=4
	#--ADD_B1=4
	#--RNA3=R1_5prime
	#--opposite-strand=TRUE
#Adapters for NEBNext Ultra II directional library prep
	#
	#
#For PE RNA-seq processed with NEBNext Ultra II Directional RNA library prep kit,
	#Read 1 contains the P5 sequence, which is the 5' end

#--ADAPT1=AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
#--ADAPT2=CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT

