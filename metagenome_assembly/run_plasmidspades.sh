#$ -S /bin/bash
#$ -N plasmidspades
#$ -V
#$ -o /workdir/users/acv46/mgmAssembly/log/pspades_$JOB_ID.out
#$ -e /workdir/users/acv46/mgmAssembly/log/pspades_$JOB_ID.err
#$ -wd /workdir/users/acv46/mgmAssembly
#$ -l h_vmem=100G
#$ -q long.q@cbsubrito
#$ -t 1-6

# last edited 13 Mar 2021, Albert Vill

WRK=/workdir/users/acv46/mgmAssembly
FQ=/home/britolab/data/PROseq/Stool
startpath=$PATH

#################################
## Switches and Global Options #############################################
#################################                                          #
                                                                           #
# increase thread restriction of cluster				   #
## default = 1								   #
## this value is used by metaSPAdes, BWA, CONCOCT, metaBAT,		   #
### MaxBin, DASTool, checkM, and GTDB-Tk 				   #
export OMP_NUM_THREADS=8						   #
									   #
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

# Create design file of file names
DESIGN_FILE=$WRK/names.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

echo -e "${NAME}: STARTING PIPELINE \n--> run < watch tail pspades_$JOB_ID.out > to view progress"

# make a directory to store all sample-specific files

PSP=$WRK/${NAME}_CAB
mkdir -p $PSP

# check if cleaned reads already exist
## if the pipeline failed post-cleaning, this saves time
## by starting with existing clean reads
## just run script again

####################
## Cleaning Reads ##
####################

echo "${NAME}: checking for existing clean reads at $FQ/${NAME}"

OUT=$PSP/fastq
mkdir -p $OUT

# edit $CLEAN1 and $CLEAN2 to tell shell where to store cleaned fastqs
# HARD-CODED DIR
CLEAN1=$FQ/${NAME}/mgm/fastq/${NAME}.clean_1.fastq
CLEAN2=$FQ/${NAME}/mgm/fastq/${NAME}.clean_2.fastq
BMTAG1=$OUT/${NAME}.clean_1.fastq
BMTAG2=$OUT/${NAME}.clean_2.fastq

if [ ! -f $BMTAG1 ] || [ ! -f $BMTAG2 ]; then

	echo "${NAME}: no existing clean reads found"
        echo "${NAME}: reading in raw fastqs"

	# unzip raw fastqs
	# note _R1/_R2 file structure
	READ1=$FQ/${NAME}_R1.fastq
	READ2=$FQ/${NAME}_R2.fastq
	gunzip ${READ1}.gz
	gunzip ${READ2}.gz

	####################
	## 1. dereplicate ##
	####################

	DEREP1=$OUT/${NAME}.derep_1.fastq
        DEREP2=$OUT/${NAME}.derep_2.fastq

	echo "${NAME}: checking for existing dereplicated reads"

	if [ ! -f $DEREP1 ] || [ ! -f $DEREP2 ]; then

		cd $OUT
		echo "${NAME}: dereplication start"

		PRINSEQ=/programs/prinseq-lite-0.20.4

		perl $PRINSEQ/prinseq-lite.pl \
			-fastq $READ1 \
			-fastq2 $READ2 \
			-derep 12345 \
			-out_format 3 \
			-out_good $OUT/${NAME}.derep \
			-out_bad $OUT/${NAME}.derep_bad

		echo "${NAME}: dereplication complete"

	else

		echo "${NAME}: dereplicated reads found, checking for trimmed reads"

	fi

	if [ ! -f $DEREP1 ] || [ ! -f $DEREP2 ]; then
    		echo "${NAME}: ERROR - dereplication failed, Aborting"
		gzip ${READ1}
        	gzip ${READ2}
		qsub -N derep_FAIL_${NAME} -M ${email} ${jobtracker}
    		exit 1
	fi

	#############
	## 2. trim ##
	#############

	TRIM1=$OUT/${NAME}.adapter_1.fastq
        TRIM2=$OUT/${NAME}.adapter_2.fastq

        if [ ! -f $TRIM1 ] || [ ! -f $TRIM2 ]; then

		echo "${NAME}: trimming start"

		TRIMMO=/programs/trimmomatic/trimmomatic-0.36.jar
		ADAPTER=$WRK/index/adapters.fasta

		java -Xmx8g -jar $TRIMMO PE \
			$DEREP1 \
			$DEREP2 \
			$OUT/${NAME}.adapter_1.fastq \
			$OUT/${NAME}.adapter_1.solo.fastq \
			$OUT/${NAME}.adapter_2.fastq \
			$OUT/${NAME}.adapter_2.solo.fastq \
			ILLUMINACLIP:${ADAPTER}:2:30:10:8:true \
			SLIDINGWINDOW:4:15 \
			LEADING:3 \
			TRAILING:3 \
			MINLEN:50

		echo "${NAME}: trimming complete"

	else

		echo "${NAME}: trimmed reads found, checking for human-depleted reads"

	fi

	if [ ! -f $TRIM1 ] || [ ! -f $TRIM2 ]; then
    		echo "${NAME}: ERROR - trimming failed, Aborting"
    		gzip ${READ1}
        	gzip ${READ2}
		qsub -N trim_FAIL_${NAME} -M ${email} ${jobtracker}
		exit 1
	fi

	###########################
	## 3. remove human reads ##
	###########################

        if [ ! -f $BMTAG1 ] || [ ! -f $BMTAG2 ]; then

		echo "${NAME}: bmtagger start"

		BMTAGGER=/programs/bmtools/bmtagger
		REFGENOME=/home/britolab/refdbs/HUMAN/Homo_sapiens_assembly19.fasta
		CONFIG=$WRK/index/bmtagger.conf
		export PATH=/programs/ncbi-blast-2.9.0+/bin:$PATH

		${BMTAGGER}/bmtagger.sh \
			-C $CONFIG \
			-b ${REFGENOME}.bitmask \
			-x ${REFGENOME}.srprism \
			-T $OUT -q1 \
			-1 $TRIM1 \
			-2 $TRIM2 \
			-o $OUT/${NAME}.clean \
			-X

		echo "${NAME}: bmtagger complete"
		echo "${NAME}: read-cleaning done"

	else

		echo "${NAME}: human-depleted reads found"
		echo "${NAME}: read-cleaning done"

	fi

	if [ ! -f $BMTAG1 ] || [ ! -f $BMTAG2 ]; then
    		echo "${NAME}: ERROR - bmtagger failed, Aborting"
		gzip ${READ1}
		gzip ${READ2}
		qsub -N bmtagger_FAIL_${NAME} -M ${email} ${jobtracker}
    		exit 1
	fi

	####################
	## 4. count reads ##
	####################

	cd $OUT

	raw1=$(echo $(wc -l ${READ1} | awk '{print $1}') / 4 | bc)
	raw2=$(echo $(wc -l ${READ2} | awk '{print $1}') / 4 | bc)
	der1=$(echo $(wc -l ${DEREP1} | awk '{print $1}') / 4 | bc)
	der2=$(echo $(wc -l ${DEREP2} | awk '{print $1}') / 4 | bc)
	trm1=$(echo $(wc -l ${TRIM1} | awk '{print $1}') / 4 | bc)
	trm2=$(echo $(wc -l ${TRIM2} | awk '{print $1}') / 4 | bc)
	bmt1=$(echo $(wc -l ${BMTAG1} | awk '{print $1}') / 4 | bc)
	bmt2=$(echo $(wc -l ${BMTAG2} | awk '{print $1}') / 4 | bc)

	# append counts from all steps
	echo -e "${NAME}\t${raw1}\t${raw2}\t${der1}\t${der2}\t${trm1}\t${trm2}\t${bmt1}\t${bmt2}" | \
	sed '1i name\traw1\traw2\tderep1\tderep2\ttrim1\ttrim2\tbmtag1\tbmtag2' \
        	> $OUT/readcounts_${NAME}_${JOB_ID}.txt

	echo -e "${NAME}: generated read counts file \n--> readcounts_${NAME}_${JOB_ID}.txt"

	##################################
	## 5. remove intermediate files ##
	##################################

	gzip ${READ1} &
	gzip ${READ2} &
	rm ${DEREP1}
	rm ${DEREP2}
	rm ${TRIM1}
	rm ${TRIM2}
	rm $OUT/${NAME}.derep_bad_1.fastq
	rm $OUT/${NAME}.derep_bad_2.fastq
	rm $OUT/${NAME}.adapter_1.solo.fastq
	rm $OUT/${NAME}.adapter_2.solo.fastq
	rm $OUT/${NAME}.derep_1_singletons.fastq
	rm $OUT/${NAME}.derep_2_singletons.fastq

	# HARD-CODED DIR
	cp ${BMTAG1} $FQ/../clean
	cp ${BMTAG2} $FQ/../clean
	cp $OUT/readcounts_${NAME}_${JOB_ID}.txt $FQ/../clean

	gzip $CLEAN1
	gzip $CLEAN2

	echo "${NAME}: checking for an existing plasmidSPAdes assembly"

else

	echo "${NAME}: clean reads found, checking for an existing plasmidSPAdes assembly"

fi

#########################
## metaSPAdes assembly ##
#########################

MSPA=$PSP/plasmidSPAdes

if [ ! -f $MSPA/contigs.fasta ]; then

	echo "${NAME}: no existing metaSPAdes assembly found"

	# if metaSPAdes previously failed, create a clean directory
	rm -r $MSPA
	mkdir -p $MSPA
	cd $MSPA

	echo "${NAME}: metaSPAdes assembly start"

	SPADES=/programs/SPAdes-3.14.0-Linux/bin
	export OMP_NUM_THREADS=8

	$SPADES/spades.py --plasmid \
		-1 $BMTAG1 \
		-2 $BMTAG2 \
		-o $MSPA \
		-t $OMP_NUM_THREADS

	if [ ! -f $MSPA/contigs.fasta ]; then
		echo "${NAME}: ERROR - metaSPAdes assembly failed, Aborting"
    		qsub -N metaSPAdes_FAIL_${NAME} -M ${email} ${jobtracker}
		exit 1
	fi

	echo "${NAME}: metaSPAdes assembly complete, checking for existing BWA alignment"

else

	echo "${NAME}: existing metaSPAdes assembly found, checking for existing BWA alignment"

fi

###################
## BWA alignment ##
###################

BWA=$MSPA/alignment
BAM=$BWA/${NAME}.bam

if [ ! -f $BAM ] || [ ! -f $BWA/${NAME}_depth.txt ]; then

	echo "${NAME}: no existing BWA alignment found"

	# if BWA previously failed, create a clean directory
        rm -r $BWA
        mkdir -p $BWA
        cd $BWA

	# create and check bwa index
	echo "${NAME}: creating BWA index of metaSPAdes assembly"

	export PATH=/programs/bwa-0.7.17:$PATH
	export PATH=/programs/samtools-1.11/bin:$PATH

	bwa index -a bwtsw -p ${NAME} $MSPA/contigs.fasta

	if [ ! -f $BWA/${NAME}.amb ] || \
	   [ ! -f $BWA/${NAME}.ann ] || \
	   [ ! -f $BWA/${NAME}.pac ] || \
	   [ ! -f $BWA/${NAME}.bwt ] || \
	   [ ! -f $BWA/${NAME}.sa ]; then
		echo "${NAME}: ERROR - BWA indexing failed, Aborting"
    		qsub -N BWAindex_FAIL_${NAME} -M ${email} ${jobtracker}
		exit 1
	fi

	# create and check bam and bam index
	echo "${NAME}: aligning reads to metaSPAdes assembly"

	bwa mem -a -v 2 -t $OMP_NUM_THREADS $BWA/${NAME} $BMTAG1 $BMTAG2 \
		> $BWA/${NAME}.unsort.bam
	samtools view -u $BWA/${NAME}.unsort.bam | \
		samtools sort -@ 4 -o $BAM

	rm $BWA/${NAME}.unsort.bam

	samtools index -@ 4 -b $BAM

	if [ ! -f $BAM ] || [ ! -f ${BAM}.bai ]; then
   		echo "${NAME}: ERROR - BWA alignment failed, Aborting"
    		qsub -N BWAmem_FAIL_${NAME} -M ${email} ${jobtracker}
		exit 1
	fi

	echo "${NAME}: BWA alignment complete"

	# create and check depth file
	DEPTH=/programs/metabat/jgi_summarize_bam_contig_depths
	echo "${NAME}: creating depth file for BWA alignment"
	$DEPTH --outputDepth $BWA/${NAME}_depth.txt $BAM

	if [ ! -f $BWA/${NAME}_depth.txt ]; then
    		echo "${NAME}: ERROR - depth calculation failed for BWA alignment, Aborting"
    		qsub -N depth_FAIL_${NAME} -M ${email} ${jobtracker}
		exit 1
	fi

	echo "${NAME}: depth calculation complete, cleaning up"

else

	echo "${NAME}: existing BWA alignment found, cleaning up"

fi

##############
## PlasFlow ##
##############

if [ ! -f $MSPA/plasflow_contigs.fasta ] || [ ! -f $MSPA/plasflow_out.tsv ]; then

	echo "${NAME}: running PlasFlow"

	export PATH=/programs/Anaconda2/bin:$PATH
	export LD_LIBRARY_PATH=/programs/Anaconda2/lib:$LD_LIBRARY_PATH
	source activate plasflow

	/programs/PlasFlow/scripts/filter_sequences_by_length.pl \
		-input $MSPA/contigs.fasta \
		-output $MSPA/plasflow_contigs.fasta \
		-thresh 1000

	PlasFlow.py \
		--input $MSPA/plasflow_contigs.fasta \
		--output $MSPA/plasflow_out.tsv \
		--threshold 0.7

	source deactivate

	if [ ! -f $MSPA/plasflow_out.tsv ]; then
                echo "${NAME}: ERROR - PlasFlow failed"
                qsub -N plasflow_FAIL_${NAME} -M ${email} ${jobtracker}
                exit 1
        fi

fi

##############
## Clean Up ##
##############

echo "${NAME} -- PIPELINE FINISHED"
echo "${NAME} -- removing intermediate bins and assembly files to save space"

qsub -N plasmids_DONE_${NAME} -M ${email} ${jobtracker}


