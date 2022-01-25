#$ -S /bin/bash
#$ -N mapHiC_redo
#$ -V
#$ -o /workdir/users/acv46/mgmAssembly/log/mapHiC_redo_$JOB_ID.out
#$ -e /workdir/users/acv46/mgmAssembly/log/mapHiC_redo_$JOB_ID.err
#$ -wd /workdir/users/acv46/mgmAssembly
#$ -l h_vmem=50G
#$ -q long.q@cbsubrito
#$ -t 1-5

# last edited 29 Mar 2021, Albert Vill

WRK=/workdir/users/acv46/mgmAssembly

#################################
## Switches and Global Options #############################################
#################################                                          #
									   #
# increase thread restriction of cluster				   #
## default = 1								   #
## this value is used by metaSPAdes, BWA, CONCOCT, metaBAT,		   #
### MaxBin, DASTool, checkM, and GTDB-Tk 				   #
export OMP_NUM_THREADS=8						   #

# to receive emails when the script fails or finishes,			   #
# add email adress here							   #
# otherwise, leave as email=NULL					   #
## NOTE -- sends one email per task					   #
email=acv46@cornell.edu							   #
									   #
# location of jobtracker.sh script					   #
jobtracker=/workdir/users/acv46/mgmAssembly/scripts/jobtracker.sh	   #
									   #
############################################################################

DESIGN_FILE=$WRK/names.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

SAMP=$WRK/${NAME}_CAB
mkdir -p $SAMP

FASTQ1_mgm=$SAMP/fastq/${NAME}.clean_1.fastq
FASTQ2_mgm=$SAMP/fastq/${NAME}.clean_2.fastq
FASTQ1_D=$SAMP/D_proxlig/${NAME}_D.clean_1.fastq
FASTQ2_D=$SAMP/D_proxlig/${NAME}_D.clean_2.fastq
FASTQ1_B=$SAMP/B_proxlig/${NAME}_B.clean_1.fastq
FASTQ2_B=$SAMP/B_proxlig/${NAME}_B.clean_2.fastq
FASTQ1_N=$SAMP/N_proxlig/${NAME}_N.clean_1.fastq
FASTQ2_N=$SAMP/N_proxlig/${NAME}_N.clean_2.fastq

echo "${NAME} -- mapping finished, N"

#####

echo "${NAME} -- starting metaphlan"

MPHL=$SAMP/metaphlan
mkdir -p $MPHL
cd $MPHL

export PATH=/programs/MetaPhlAn-2.0:/programs/MetaPhlAn-2.0/utils:$PATH

metaphlan2.py \
        --input_type fastq \
        --bowtie2db /workdir/users/fnn3/references/db_v20/mpa_v20_m200 \
        --bowtie2out ${NAME}_D.bowtie2.bz2 \
        --nproc $OMP_NUM_THREADS \
        $FASTQ1_D,$FASTQ2_D \
        ${NAME}_D_profile.txt

metaphlan2.py \
        --input_type fastq \
        --bowtie2db /workdir/users/fnn3/references/db_v20/mpa_v20_m200 \
        --bowtie2out ${NAME}_B.bowtie2.bz2 \
        --nproc $OMP_NUM_THREADS \
        $FASTQ1_B,$FASTQ2_B \
        ${NAME}_B_profile.txt

metaphlan2.py \
        --input_type fastq \
        --bowtie2db /workdir/users/fnn3/references/db_v20/mpa_v20_m200 \
        --bowtie2out ${NAME}_N.bowtie2.bz2 \
        --nproc $OMP_NUM_THREADS \
        $FASTQ1_N,$FASTQ2_N \
        ${NAME}_N_profile.txt
