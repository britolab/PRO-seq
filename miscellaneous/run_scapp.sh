#$ -S /bin/bash
#$ -N scapp
#$ -V
#$ -o /workdir/users/acv46/mgmAssembly/log/scapp_$JOB_ID.out
#$ -e /workdir/users/acv46/mgmAssembly/log/scapp_$JOB_ID.err
#$ -wd /workdir/users/acv46/mgmAssembly
#$ -l h_vmem=50G
#$ -t 1-2
#$ -q long.q@cbsubrito

# uses existing metaSPAdes graph and bam alignment from CAB pipeline

WRK=/workdir/users/acv46/mgmAssembly
export OMP_NUM_THREADS=16

# see /home/britolab/acv46/SCAPP/_README for install details
source /home/britolab/acv46/miniconda3/bin/activate
conda activate scapp

list=$WRK/samples.txt
List=${list}
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $List)
        samp=`basename "$DESIGN"`

OUT=$WRK/${samp}_CAB/scapp
mkdir -p $OUT

graph=$WRK/${samp}_CAB/metaSPAdes/assembly_graph.fastg
bam=$WRK/${samp}_CAB/metaSPAdes/alignment/${samp}.bam

echo "${samp} -- starting scapp"

scapp \
	-g $graph \
	-o $OUT \
	-b $bam \
	-p $OMP_NUM_THREADS \
	-k 55 \
	-l 1000 \
	-pst 0.9

echo "${samp} -- scapp finished"

conda deactivate
