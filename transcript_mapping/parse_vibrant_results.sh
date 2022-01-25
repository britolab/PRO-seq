#$ -S /bin/bash
#$ -N parse_vibrant
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq2/log/parse_vibrant_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq2/log/parse_vibrant_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=10G
#$ -q long.q@cbsubrito
#$ -t 1-2

WRK=/workdir/users/acv46/stool_PROSeq2/phage/vibrant
samples=$WRK/../samples.txt
CDHIT=/programs/cd-hit-4.8.1/cd-hit-est

DESIGN_FILE=${samples}
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

base=$(echo ${NAME} | awk -F_ '{print $1}')

cd $WRK/${base}

lyt_ffn=VIBRANT_*/VIBRANT_phages_*/*.phages_lytic.ffn
lys_ffn=VIBRANT_*/VIBRANT_phages_*/*.phages_lysogenic.ffn

# terminase large subunit

grep -A1 -i --no-filename "terminase" $lyt_ffn | \
	grep -A1 -i "large" | \
	grep -v -- "^--$" \
	> lytic_terminase_large.fa

grep -A1 -i --no-filename "terminase" $lys_ffn | \
        grep -A1 -i "large" | \
        grep -v -- "^--$" \
        > lysogenic_terminase_large.fa

# terminase small subunit

grep -A1 -i --no-filename "terminase" $lyt_ffn | \
        grep -A1 -i "small" | \
        grep -v -- "^--$" \
        > lytic_terminase_small.fa

grep -A1 -i --no-filename "terminase" $lys_ffn | \
        grep -A1 -i "large" | \
        grep -v -- "^--$" \
        > lysogenic_terminase_small.fa

# major capsid proteins

grep -A1 -i --no-filename "capsid" $lyt_ffn | \
        grep -v -- "^--$" \
        > lytic_capsid.fa

grep -A1 -i --no-filename "capsid" $lys_ffn | \
        grep -v -- "^--$" \
        > lysogenic_capsid.fa

# tail fibers

grep -A1 -i --no-filename "tail fiber" $lyt_ffn | \
        grep -v -- "^--$" \
        > lytic_tailfiber.fa

grep -A1 -i --no-filename "tail fiber" $lys_ffn | \
        grep -v -- "^--$" \
        > lysogenic_tailfiber.fa

# CI repressors

grep -A1 -i --no-filename " ci" $lyt_ffn | \
        grep -vi "circularization" | \
        grep -v "^--$" \
        > lytic_ci_repressor.fa

grep -A1 -i --no-filename " ci" $lys_ffn | \
        grep -vi "circularization" | \
        grep -v "^--$" \
        > lysogenic_ci_repressor.fa
