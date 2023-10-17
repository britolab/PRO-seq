#$ -S /bin/bash
#$ -N proseq2.0
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq3/log/proseq_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq3/log/proseq_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq3
#$ -l h_vmem=70G
#$ -q long.q@cbsubrito
#$ -t 1-22
#$ -tc 4
#$ -pe parenv 8

# corrected orientation based on observed tRNA mapping
# --opposite-strand=TRUE doesn't do anything

# this is the caller script for the Danko Lab script (proseq2.0_revised.bsh)
# only map to contigs > 1000 bp

# stackoverflow.com/questions/10496758

# PRO-seq script dependencies all installed in one conda env


WRK=/workdir/users/acv46/stool_PROSeq3
OUT=$WRK/transcriptomes
FQ=/home/britolab/data/NextSeq_10458816

mkdir -p $OUT

#generate bwaIndex using "bwa index" (generated automatically during metaSPADES assembly)
#generate chromInfo file using Get_Lengths.sh

export PATH=/programs/kentUtils/bin:$PATH
export PATH=/programs/seqtk:$PATH
export PATH=/programs/prinseq-lite-0.20.4:$PATH
export PATH=/programs/proseq2.0:$PATH
export PATH=/programs/cutadapt-4.1/bin:$PATH
export PYTHONPATH=/programs/cutadapt-4.1/lib/python3.9/site-packages:/programs/cutadapt-4.1/lib64/python3.9/site-packages

export PATH=/programs/bwa-0.7.17:$PATH
export PATH=/programs/samtools-1.15/bin:$PATH
export PATH=/programs/bedtools2-2.29.2/bin:$PATH
export PATH=/programs/bedops-2.4.35/bin:$PATH
export OMP_NUM_THREADS=8

LIST=$OUT/names.txt
  DESIGN=$(sed -n "${SGE_TASK_ID}p" $LIST)
  NAME=`basename "$DESIGN"`

raw1=$(grep -P "\t${NAME}" $OUT/names_to_reads.txt | cut -f1 | head -n1)
raw2=$(grep -P "\t${NAME}" $OUT/names_to_reads.txt | cut -f1 | tail -n1)

read1=${FQ}/${raw1}
read2=${FQ}/${raw2}

echo "${NAME} -- starting proseq2.0 pipeline"
echo "        --> read1 = ${read1}"
echo "        --> read2 = ${read2}"

# exit for debugging
# exit 1

base=$(echo ${NAME} | cut -f1 -d'_')
MGM=$WRK/assembly/${base}_25May2022
samp=$OUT/${base}
mkdir -p ${samp}/index

cp ${read1} ${samp}
mv ${samp}/${raw1} ${samp}/${NAME}_R1.fastq.gz
cp ${read2} ${samp}
mv ${samp}/${raw2} ${samp}/${NAME}_R2.fastq.gz

# subset contigs by length

len=1000
allcont=$MGM/metaSPAdes/contigs.fasta
subcont=${samp}/index/contigs_${len}.fasta

if [ ! -f ${subcont} ]; then
  echo "${NAME} -- subsetting contigs to length ${len}"
  awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' ${allcont} | \
    awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' | \
    awk -v len=${len} -F$'\t' '$1 >= len' | \
    sort -k1,1nr | \
    cut -f 2- | \
    tr "\t" "\n" \
    > ${subcont}
else
  echo "${NAME} -- subset contigs found"
fi

# access bwa index

echo "${NAME} -- checking for bwa index"

bwaIndex=${samp}/index/${base}

if [ ! -f ${bwaIndex}.amb ] || \
   [ ! -f ${bwaIndex}.ann ] || \
   [ ! -f ${bwaIndex}.pac ] || \
   [ ! -f ${bwaIndex}.bwt ] || \
   [ ! -f ${bwaIndex}.sa ]; then
  echo "${NAME}: creating bwa index"
  cd ${samp}/index
  bwa index -a bwtsw -p ${base} ${subcont}
else
  echo "${NAME}: bwa index found"
fi

# create chromInfo file

cinfo=${samp}/index/chromInfo.txt
echo "${NAME} -- checking for chromInfo file (min length ${len} bp)"

if [ ! -f ${cinfo} ]; then
  awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' ${subcont} > ${cinfo}
else
  echo "${NAME}: chromInfo.txt found"
fi

# process and map reads

TMP=${samp}/temp
mkdir -p $TMP

if [[ "$NAME" == *"PRO"* ]]; then

  echo "${NAME} -- mapping PRO-seq reads"

  cd ${samp}

  proseq2.0.bsh \
    -i $bwaIndex \
    -c $cinfo \
    -PE \
    --RNA3=R2_5prime \
    --RNA5=R1_5prime \
    --ADAPT2=GATCGTCGGACTGTAGAACTCTGAAC \
    --ADAPT1=TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC \
    --UMI1=0 \
    --UMI2=6 \
    --ADD_B1=6 \
    --ADD_B2=6 \
    -T $TMP \
    -O ${samp} \
    -I ${NAME} \
    --thread=${OMP_NUM_TREADS}

  # bash $PRO \
  #   -i $bwaIndex \
  #   -c $cinfo \
  #   -PE \
  #   --RNA3=R1_5prime \
  #   --RNA5=R2_5prime \
  #   --ADAPT1=GATCGTCGGACTGTAGAACTCTGAAC \
  #   --ADAPT2=TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC \
  #   --UMI1=6 \
  #   --UMI2=0 \
  #   --ADD_B1=6 \
  #   --ADD_B2=6 \
  #   -T $TMP \
  #   -O ${samp} \
  #   -I ${NAME} \
  #   --thread=${OMP_NUM_TREADS} \
  #   --opposite-strand=TRUE

elif [[ "$NAME" == *"RNA"* ]] || \
     [[ "$NAME" == *"TEX"* ]]; then

  echo "${NAME} -- mapping NEBNext-processed reads"

  cd ${samp}

  proseq2.0.bsh \
    -i $bwaIndex \
    -c $cinfo \
    -PE \
    --RNA3=R2_5prime \
    --RNA5=R1_5prime \
    --ADAPT2=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --ADAPT1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --UMI2=11 \
    --UMI1=0 \
    --ADD_B1=8 \
    --ADD_B2=8 \
    -T $TMP \
    -O ${samp} \
    -I ${NAME} \
    --thread=${OMP_NUM_TREADS}

  # bash $PRO \
  #   -i $bwaIndex \
  #   -c $cinfo \
  #   -PE \
  #   --RNA3=R1_5prime \
  #   --RNA5=R2_5prime \
  #   --ADAPT1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  #   --ADAPT2=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  #   --UMI1=11 \
  #   --UMI2=0 \
  #   --ADD_B1=8 \
  #   --ADD_B2=8 \
  #   -T $TMP \
  #   -O ${samp} \
  #   -I ${NAME} \
  #   --thread=${OMP_NUM_TREADS} \
  #   --opposite-strand=TRUE

fi

echo "${NAME} -- removing copied fastqs"

rm ${samp}/${NAME}_R1.fastq.gz
rm ${samp}/${NAME}_R2.fastq.gz

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

