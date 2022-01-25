#$ -S /bin/bash
#$ -N vs2
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq2/log/vs2_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq2/log/vs2_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=100G
#$ -q long.q@cbsubrito
#$ -t 1-2

# installation
# conda install virsorter
# conda update virsorter
# virsorter setup -d db -j 4


WRK=/workdir/users/acv46/stool_PROSeq2/phage/virsorter
VSORT=/home/acv46/miniconda3/bin/virsorter
export OMP_NUM_THREADS=16

DESIGN_FILE=$WRK/assemblies.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

echo "starting ${NAME}"

source /home/acv46/miniconda3/bin/activate

fasta=/workdir/users/acv46/mgmAssembly/${NAME}_CAB/metaSPAdes/contigs.fasta

cd $WRK

cp $fasta $WRK/${NAME}.fasta

$VSORT run -w ${NAME}.out -i ${NAME}.fasta -j 16



