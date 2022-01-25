#$ -S /bin/bash
#$ -N gtdbtk
#$ -V
#$ -o /workdir/users/acv46/mgmAssembly/log/gtdbtk_$JOB_ID.out
#$ -e /workdir/users/acv46/mgmAssembly/log/gtdbtk_$JOB_ID.err
#$ -wd /workdir/users/acv46/mgmAssembly
#$ -l h_vmem=160G
#$ -q long.q@cbsubrito
#$ -t 1

# This script runs vamb with different hyperparameter settings
# uses existing metaSPAdes alignments from DASTool pipeline

WRK=/workdir/users/acv46/mgmAssembly/vamb/out_465249
MDM=$WRK/medium_NN
SML=$WRK/small_NN
LRG=$WRK/large_NN

export PYTHONPATH=/programs/gtdbtk-1.0.2/lib/python3.6/site-packages/
export PATH=/programs/gtdbtk-1.0.2/bin:$PATH
export PATH=/programs/hmmer/binaries:$PATH
export PATH=/programs/prodigal-2.6.3:$PATH
export PATH=/programs/FastTree-2.1.10:$PATH
export PATH=/programs/fastANI-1.3:$PATH
export PATH=/programs/pplacer-Linux-v1.1.alpha19:$PATH
export GTDBTK_DATA_PATH=/workdir/GtdbTK/release89/

mkdir -p $MDM/gtdbtk

gtdbtk classify_wf \
	--cpus 8 \
        --genome_dir $MDM/bins \
        --extension fna \
        --out_dir $MDM/gtdbtk

mkdir -p $SML/gtdbtk

gtdbtk classify_wf \
        --cpus 8 \
        --genome_dir $SML/bins \
        --extension fna \
        --out_dir $SML/gtdbtk

mkdir -p $LRG/gtdbtk

gtdbtk classify_wf \
        --cpus 8 \
        --genome_dir $LRG/bins \
        --extension fna \
        --out_dir $LRG/gtdbtk
