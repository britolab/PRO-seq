#$ -S /bin/bash
#$ -N parse_vibrant
#$ -V
#$ -o /workdir/users/acv46/ben_May2021/log/parse_vibrant_ci_$JOB_ID.out
#$ -e /workdir/users/acv46/ben_May2021/log/parse_vibrant_ci_$JOB_ID.err
#$ -wd /workdir/users/acv46/ben_May2021
#$ -l h_vmem=10G
#$ -q long.q@cbsubrito2
#$ -t 1

WRK=/workdir/users/acv46/ben_May2021/phage/vibrant
samples=/workdir/users/acv46/ben_May2021/vamb/small_NN/samples.txt
export OMP_NUM_THREADS=8
CDHIT=/programs/cd-hit-4.8.1/cd-hit-est

#DESIGN_FILE=${samples}
#        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
#        NAME=`basename ${DESIGN}`

cd $WRK

lyt_faa=$WRK/*/VIBRANT_*_merged/VIBRANT_phages_*_merged/*_merged.phages_lytic.faa
lys_faa=$WRK/*/VIBRANT_*_merged/VIBRANT_phages_*_merged/*_merged.phages_lysogenic.faa
lyt_ffn=$WRK/*/VIBRANT_*_merged/VIBRANT_phages_*_merged/*_merged.phages_lytic.ffn
lys_ffn=$WRK/*/VIBRANT_*_merged/VIBRANT_phages_*_merged/*_merged.phages_lysogenic.ffn

# tail fibers

grep -A1 -i --no-filename " ci" $lyt_ffn | \
	grep -vi "circularization" | \
        grep -v "^--$" \
        > lytic_ci_repressor.fa

grep -A1 -i --no-filename " ci" $lys_ffn | \
        grep -vi "circularization" | \
	grep -v "^--$" \
        > lysogenic_ci_repressor.fa
