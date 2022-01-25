#$ -S /bin/bash
#$ -N vamb
#$ -V
#$ -o /workdir/users/acv46/mgmAssembly/log/vamb_$JOB_ID.out
#$ -e /workdir/users/acv46/mgmAssembly/log/vamb_$JOB_ID.err
#$ -wd /workdir/users/acv46/mgmAssembly
#$ -l h_vmem=100G
#$ -q long.q@cbsubrito
#$ -t 1

# This script runs vamb with different hyperparameter settings
# uses existing metaSPAdes alignments from DASTool pipeline

SPA=/workdir/users/acv46/mgmAssembly/metaSPAdes/*/contigs.fasta
MAP=/programs/minimap2-2.17/minimap2
FQ=/workdir/users/acv46/mgmAssembly/fastq
WRK=/workdir/users/acv46/mgmAssembly/vamb
OUT=$WRK/out_${JOB_ID}
BAM=$OUT/bam
MDM=$OUT/medium_NN
SML=$OUT/small_NN
LRG=$OUT/large_NN

mkdir -p $OUT

export PATH=$HOME/.local/bin:$PATH
export PYTHONPATH=$HOME/.local/lib/python3.9/site-packages

# uses existing metaSPAdes alignments from DASTool pipeline

# concatenate assemblies into single file
## use -m parameter to discard small sequences (default 2000)

echo "${JOB_ID}: running concatenate.py"

concatenate.py $OUT/catalog.fna.gz $SPA

# minimap2 alignment
## Note: vamb requires bam sorted by name
## However, bwa cannot create an index for anything other than\
## coordinate-sorted bam files
## just use minimap2

# make index of contig catalog

echo "${JOB_ID}: creating minimap2 index of contig catalog"

$MAP -d $OUT/catalog.mmi $OUT/catalog.fna.gz

# get sample names from $FQ
## assumes [SAMPLE.clean_1.fastq] and [SAMPLE.clean_2.fastq] file structure
## assumes all files in $FQ correspond to assemblies in $SPA

find $FQ -type f -iname "*.clean_1.fastq" | \
awk -F "/" '{print $NF}' | \
sed 's/.clean_1.fastq//g' > $OUT/samples.txt

## lots of options
## see $MAP --help

echo "${JOB_ID}: mapping $(wc -l < $OUT/samples.txt) samples to contig catalog"

mkdir -p $BAM

while read line; do

	echo "${JOB_ID}: mapping sample $(grep -n "${line}" $OUT/samples.txt | awk -F ":" '{print $1}') of $(wc -l < $OUT/samples.txt) -- ${line}"

	$MAP \
	-t 16 \
	-N 50 \
	-ax sr \
	$OUT/catalog.mmi \
	$FQ/${line}.clean_1.fastq \
	$FQ/${line}.clean_2.fastq | \
	samtools view \
	-F 3584 \
	-b \
	--threads 8 > $BAM/${line}.bam

done < $OUT/samples.txt

# vamb

echo "${JOB_ID}: starting vamb with default hyperparameters"
## Default hyperparameters
vamb \
-l 32 \
-n 512 512 \
--outdir $MDM \
--fasta $OUT/catalog.fna.gz \
--bamfiles $BAM/*.bam \
-o C \
--minfasta 100000

echo "${JOB_ID}: starting vamb with decreased hyperparameters"
## Smaller Network
vamb \
-l 24 \
-n 384 384 \
--outdir $SML \
--fasta $OUT/catalog.fna.gz \
--bamfiles $BAM/*.bam \
-o C \
--minfasta 100000

echo "${JOB_ID}: starting vamb with increased hyperparameters"
## Larger Network
vamb \
-l 40 \
-n 768 768 \
--outdir $LRG \
--fasta $OUT/catalog.fna.gz \
--bamfiles $BAM/*.bam \
-o C \
--minfasta 100000


#### CheckM prerequsites and installation

# Install software in your home directory($HOME/.local/). (You only need to do it once)

#pip3 install --user checkm-genome

# Set environment (You need to run it everytime you create a new session)
export PATH=$HOME/.local/bin:$PATH
export PYTHONPATH=$HOME/.local/lib/python3.6/site-packages
export PATH=/programs/hmmer/bin:$PATH
export PATH=/programs/prodigal-2.6.3:$PATH
export PATH=/programs/pplacer-Linux-v1.1.alpha19:$PATH

echo "${JOB_ID}: running checkM on vamb bins -- default hyperparameters"

mkdir -p $MDM/checkm
cd $MDM/checkm

checkm lineage_wf -q \
	-t 10 \
	--pplacer_threads 16 \
	-x fna \
	-f checkM_lineage_M.txt \
	$MDM/bins \
	$MDM/checkm

echo "${JOB_ID}: running checkM on vamb bins -- decreased hyperparameters"

mkdir -p $SML/checkm
cd $SML/checkm

checkm lineage_wf -q \
        -t 10 \
        --pplacer_threads 16 \
        -x fna \
        -f checkM_lineage_S.txt \
        $SML/bins \
        $SML/checkm

echo "${JOB_ID}: running checkM on vamb bins -- increased hyperparameters"

mkdir -p $LRG/checkm
cd $LRG/checkm

checkm lineage_wf -q \
        -t 10 \
        --pplacer_threads 16 \
        -x fna \
        -f checkM_lineage_L.txt \
        $LRG/bins \
        $LRG/checkm

echo "${JOB_ID}: vamb script done"
