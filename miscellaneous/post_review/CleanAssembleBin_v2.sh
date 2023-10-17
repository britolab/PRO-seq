#$ -S /bin/bash
#$ -N CABv2
#$ -V
#$ -o /workdir/users/acv46/stool_PROSeq3/log/CABv2_$JOB_ID.out
#$ -e /workdir/users/acv46/stool_PROSeq3/log/CABv2_$JOB_ID.err
#$ -wd /workdir/users/acv46/stool_PROSeq3/
#$ -l h_vmem=250G
#$ -q long.q@cbsubrito
#$ -t 1-2
#$ -tc 1
#$ -pe parenv 16

# last edited 31 Oct 2022, Albert Vill
# this is the updated version of the CAB pipeline
# includes Rosella binning and bin summary script

##########################################
## CAB Pipeline -- Clean, Assemble, Bin ##
##########################################

# For raw paired-end Illumina reads, this script does the following:
####################################################################
## dereplicates reads with prinseq
## trims reads with Trimmomatic [consider using BBDuk instead]
## removes reads mapping to human genome with bmtagger
## counts reads at each step
## assembles reads with metaSPAdes
## aligns reads to metaSPAdes assembly to get bam
## bins contigs with CONCOCT, MetaBAT, MaxBin, and Rosella
## merges bins with DASTool
## evaluates DASTool bins with CheckM
## (optional) calls taxa on bins with gtdbtk
#### Note, gtdbtk uses pplacer, which requires a lot of memory (>150GB)
#### if you choose to run gtdbtk (switch below), increase h_vmem to 150GB
## calls operons with POEM

# Before running:
#################
## edit qsub parameters -o, -e, -wd above (out file, error file, workdir)
## designate a working directory as $WRK, below
## create a file named 'names.txt' in working dir with sample names
## designate directory containing raw reads as $FQ
#### Note, the $FQ directory structure is hard-coded as $FQ/${NAME}/mgm/fastq
#### see comments in the code (# HARD-CODED DIR) to see where directories need to be edited
## edit the -t parameter above to specify the number of jobs (e.g. for 4 samples, '-t 1-4')
## edit $READ1 and $READ2 to point to your raw fastq.gz files

# Software dependencies and locations:
######################################
## prinseq	--      /programs/prinseq-lite-0.20.4
## trimmomatic	--      /programs/trimmomatic/trimmomatic-0.36.jar
## BLAST	--      /programs/ncbi-blast-2.9.0+/bin
## bmtagger	--      /programs/bmtools/bmtagger
## SPAdes	--      /programs/SPAdes-3.14.0-Linux/bin
## BWA		--	/programs/bwa-0.7.17
## samtools	--	/programs/samtools-1.11/bin
## CONCOCT	--	/programs/miniconda3/bin/activate concoct
## metaBAT2	--	/programs/metabat
## MaxBin	--	/programs/MaxBin-2.2.4
## DASTool	--	/home/acv46/DAS_Tool
## diamond	--	/programs/diamond-2.0.4
## usearch	--	/programs/usearch11.0.667
## ruby		--	/programs/ruby/bin
## pullseq	--	/home/britolab/acv46/pullseq/bin
## R		--	/programs/R-4.0.0/bin
## CheckM	--	/home/acv46/miniconda3/envs/checkm
## CheckM2	--	/home/acv46/miniconda3/envs/checkm2
## Rosella	--	/home/acv46/miniconda3/envs/rosella
## GTDB-Tk	--	/programs/gtdbtk-2.1.1/bin
## prodigal	--	/programs/prodigal-2.6.3
## HMMER	--	/programs/hmmer-3.3/bin
## FastTree	--	/programs/FastTree-2.1.10
## pplacer	--	/programs/pplacer-Linux-v1.1.alpha19
## FastANI	--	/programs/fastANI-1.3
## POEM		--	/home/britolab/acv46/POEM_py3k

# To Do:
########
## create a summary file of all relevant values
### fastq: input, derep, trim, and bmtag reads
### bwa: N50 and number of contigs > 1000 bp
### bins: number of bins for each program
### checkM: nubmer of HQ bins
### poem: number and size of operons
## fix checkM plotting dependencies
## replace POEM with metaRON
## eliminate redundant gene annotations

# Set directories to read/write files:
######################################

# HARD-CODED DIR
## working directory where sample-specific subdirs will be made
## index folder goes here, too
WRK=/workdir/users/acv46/stool_PROSeq3/assembly
## directory containing raw fastqs
FQ=/home/britolab/data/US_mgms
startpath=$PATH

#################################
## Switches and Global Options #############################################
#################################                                          #
                                                                           #
# Switch for using GTDB-Tk						   #
## 1 = run GTDB-Tk, 0 = don't run					   #
## If 1, allocate at least 150 GB memory per qsub task (-l h_vmem=150G)    #
## note: GTDB-Tk doesn't work reliably on brito1 \                         #
## copy bins to brito2 and run there, specifying export OMP_NUM_THREADS=8  #
gFlag=1									   #
									   #
# Switch for calling operons with POEM					   #
## 1 = run POEM, 0 = don't						   #
pFlag=0									   #
									   #
# Switch for MetaRon gene prediction tool				   #
## 1 = MetaGeneMark, 2 = Prodigal					   #
mTool=2									   #
									   #
# increase thread restriction of cluster				   #
## default = 1								   #
## this value is used by metaSPAdes, BWA, CONCOCT, metaBAT,		   #
### MaxBin, DASTool, checkM, and GTDB-Tk 				   #
export OMP_NUM_THREADS=16						   #
									   #
# DASTool parameters							   #
## search engine can be diamond, blast, or usearch			   #
searcheng=diamond							   #
## score threshold							   #
score=0									   #
									   #
# Completeness and Contamination cutoffs for DASTool bins		   #
## only bins passing these checkM thresholds will be processed by POEM     #
complete=90								   #
contam=10								   #
									   #
# to receive emails when the script fails or finishes,			   #
# add email adress here							   #
# otherwise, leave as email=NULL					   #
## NOTE -- sends one email per task					   #
email=acv46@cornell.edu							   #
									   #
# location of jobtracker.sh script					   #
jobtracker=/workdir/users/acv46/stool_PROSeq3/scripts/jobtracker.sh	   #
									   #
############################################################################

# Create design file of file names
DESIGN_FILE=$WRK/names.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename ${DESIGN}`

echo -e "${NAME}: STARTING PIPELINE \n--> run < watch tail CAB_$JOB_ID.out > to view progress"

# make a directory to store all sample-specific files

CAB=$WRK/${NAME}
mkdir -p $CAB

# check if cleaned reads already exist
## if the pipeline failed post-cleaning, this saves time
## by starting with existing clean reads
## just run script again

####################
## Cleaning Reads ##
####################

echo "${NAME}: checking for existing clean reads at $FQ/${NAME}"

OUT=$CAB/fastq
mkdir -p $OUT

# edit $CLEAN1 and $CLEAN2 to tell shell where to store cleaned fastqs
# HARD-CODED DIR
CLEAN1=$FQ/${NAME}/clean/${NAME}.clean_1.fastq
CLEAN2=$FQ/${NAME}/clean/${NAME}.clean_1.fastq
BMTAG1=$OUT/${NAME}.clean_1.fastq
BMTAG2=$OUT/${NAME}.clean_2.fastq

if [ ! -f $BMTAG1 ] || [ ! -f $BMTAG2 ]; then

	echo "${NAME}: no existing clean reads found"
        echo "${NAME}: reading in raw fastqs"

	# unzip raw fastqs
	# HARD-CODED DIR
	# note _R1/_R2 file structure
	READ1=$FQ/${NAME}/raw/${NAME}_1.fq
	READ2=$FQ/${NAME}/raw/${NAME}_2.fq
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
		export PATH=/programs/ncbi-blast-2.12.0+/bin:$PATH

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
	cp ${BMTAG1} $FQ/${NAME}/clean
	cp ${BMTAG2} $FQ/${NAME}/clean
	cp $OUT/readcounts_${NAME}_${JOB_ID}.txt $FQ/${NAME}/mgm/fastq

	gzip $CLEAN1
	gzip $CLEAN2

	echo "${NAME}: checking for an existing metaSPAdes assembly"

else

	echo "${NAME}: clean reads found, checking for an existing metaSPAdes assembly"

fi

#########################
## metaSPAdes assembly ##
#########################

MSPA=$CAB/metaSPAdes

if [ ! -f $MSPA/contigs.fasta ]; then

	echo "${NAME}: no existing metaSPAdes assembly found"

	# if metaSPAdes previously failed, create a clean directory
	rm -r $MSPA
	mkdir -p $MSPA
	cd $MSPA

	echo "${NAME}: metaSPAdes assembly start"

	SPADES=/programs/SPAdes-3.15.4-Linux/bin

	# make sure -m is equal to or less than value specified by qsub -l param
	$SPADES/metaspades.py \
		-1 $BMTAG1 \
		-2 $BMTAG2 \
		-o $MSPA \
		-t $OMP_NUM_THREADS \
		-m 250

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
	export PATH=/programs/samtools-1.15.1/bin:$PATH

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
		samtools sort -@ $OMP_NUM_THREADS -o $BAM

	rm $BWA/${NAME}.unsort.bam

	samtools index -@ $OMP_NUM_THREADS -b $BAM

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

	echo "${NAME}: depth calculation complete, checking for existing DASTool bins"

else

	echo "${NAME}: existing BWA alignment found, checking for existing DASTool bins"

fi

DAS=$CAB/DASTool
CONC=$CAB/CONCOCT
BAT=$CAB/metaBAT
MXB=$CAB/MaxBin

# set seed for metaBAT
seed=123

if [ $(find $DAS/${NAME}_DASTool_bins -name '*.fa' | wc -l) -eq 0 ]; then

	echo "${NAME}: no existing DASTool bins found, starting binning pipeline"

        mkdir -p $DAS

	#####################
	## CONCOCT binning ##
	#####################

	echo "${NAME}: checking for existing CONCOCT bins"

	if [ $(find $CONC -name '*.fa' | wc -l) -eq 0 ]; then

		echo "${NAME}: no existing CONCOCT bins found"

		# if CONCOCT previously failed, create a clean directory
        	rm -r $CONC
        	mkdir -p $CONC
        	cd $CONC

		source /programs/miniconda3/bin/activate concoct

		#### Step 1: Cut contigs into smaller parts
		echo "${NAME}: CONCOCT, Cut up fasta"
		cut_up_fasta.py $MSPA/contigs.fasta \
			-c 10000 \
			-o 0 \
			--merge_last \
			-b $CONC/${NAME}_contigs.10k.bed \
			> $CONC/${NAME}_contigs.10k.fa

		#### Step 2: Generate table with coverage depth info per sample and subcontig
		echo "${NAME}: CONCOCT, Coverage table start"
		concoct_coverage_table.py $CONC/${NAME}_contigs.10k.bed $BAM \
			> $CONC/${NAME}_coverage_table.tsv

		#### Step 3: Run CONCOCT
		echo "${NAME}: CONCOCT binning start"
		concoct \
			--composition_file $CONC/${NAME}_contigs.10k.fa \
			--coverage_file $CONC/${NAME}_coverage_table.tsv \
			-b $CONC \
			-t $OMP_NUM_THREADS

		#### Step 4: Merge subcontig clustering into original contig clustering
		echo "${NAME}: CONCOCT, merge cutup clustering start"
		merge_cutup_clustering.py $CONC/clustering_gt1000.csv \
			> $CONC/clustering_merged.csv

		#### Step 5: Extract bins as individual FASTA
		echo "${NAME}: CONCOCT, extract bins start"
		extract_fasta_bins.py \
			$MSPA/contigs.fasta \
			$CONC/clustering_merged.csv \
			--output_path $CONC

		echo "${NAME}: CONCOCT finished"

		conda deactivate

		echo "${NAME}: CONCOCT detected $(find $CONC -name '*.fa' | wc -l) bins"

	else

		echo "${NAME}: found $(find $CONC -name '*.fa' | wc -l) existing CONCOCT bins"

	fi

	#####################
	## metaBAT binning ##
	#####################

	echo "${NAME}: checking for existing metaBAT bins"

	if [ $(find $BAT/contigs.fasta.metabat-bins${seed} -name '*.fa' | wc -l) -eq 0 ]; then

		echo "${NAME}: no existing metaBAT bins found"

		# if metaBAT previously failed, create a clean directory
        	rm -r $BAT
        	mkdir -p $BAT
	        cd $BAT

		echo "${NAME}: metaBAT start"
		/programs/metabat/runMetaBat.sh \
			-m 1500 \
			-s 10000 \
			-t ${OMP_NUM_THREADS} \
			--seed ${seed} \
			$MSPA/contigs.fasta ${BAM}

		echo "${NAME}: metaBAT finished"
		echo "${NAME}: metaBAT detected $(find $BAT/contigs.fasta.metabat-bins${seed} -name '*.fa' | wc -l) bins"

	else

		echo "${NAME}: found $(find $BAT/contigs.fasta.metabat-bins${seed} -name '*.fa' | wc -l) existing metaBAT bins"

	fi

	####################
	## MaxBin binning ##
	####################

	echo "${NAME}: checking for existing MaxBin bins"

	if [ $(find $MXB -name '*.fasta' | wc -l) -eq 0 ]; then

		echo "${NAME}: no existing MaxBin bins found"

		# if MaxBin previously failed, create a clean directory
        	rm -r $MXB
        	mkdir -p $MXB
        	cd $MXB

		echo "${NAME}: MaxBin start"
		/programs/MaxBin-2.2.4/run_MaxBin.pl \
			-contig $MSPA/contigs.fasta \
			-reads $BMTAG1 \
			-reads2 $BMTAG2 \
			-max_iteration 50 \
			-thread ${OMP_NUM_THREADS} \
			-out ${NAME}

		echo "${NAME}: MaxBin finished"
		echo "${NAME}: MaxBin detected $(find $MXB -name '*.fasta' | wc -l) bins"

	else

		echo "${NAME}: found $(find $MXB -name '*.fasta' | wc -l) existing MaxBin bins"

	fi

	#####################
        ## Rosella binning ##
        #####################

	

	#############
	## DASTool ##
	#############

	dastool=/home/britolab/acv46/DAS_Tool

        echo "${NAME}: DASTool, fasta to scaffolds"

        export PATH=/programs/ncbi-blast-2.12.0+/bin:$PATH
        export PATH=/programs/diamond-2.0.9:$PATH
        export PATH=/programs/usearch11.0.667:$PATH
        export PATH=/programs/prodigal-2.6.3:$PATH
        export PATH=/programs/ruby/bin:$PATH
        export PATH=/home/britolab/acv46/pullseq/bin:$PATH
	export PATH=/programs/R-4.0.0/bin:$PATH

	echo "$DAS"

        if [ ! -f $DAS/${NAME}_maxbin2_scaffolds2bin.tsv ]; then

                echo "${NAME}: running fasta to scaffolds, MaxBin"
                $dastool/src/Fasta_to_Scaffolds2Bin.sh -i $MXB -e fasta \
                        > $DAS/${NAME}_maxbin2_scaffolds2bin.tsv

        else

		echo "${NAME}: using existing file for MaxBin: $DAS/${NAME}_maxbin2_scaffolds2bin.tsv"

        fi

        # note the metabat seed number
        if [ ! -f $DAS/${NAME}_metabat2_scaffolds2bin.tsv ]; then

                echo "${NAME}: running fasta to scaffolds, metaBAT"
                $dastool/src/Fasta_to_Scaffolds2Bin.sh -i $BAT/contigs.fasta.metabat-bins${seed} -e fa \
                        > $DAS/${NAME}_metabat2_scaffolds2bin.tsv

        else

                echo "${NAME}: using existing file for metaBAT: $DAS/${NAME}_metabat2_scaffolds2bin.tsv"

        fi

	if [ ! -f $DAS/${NAME}_concoct_scaffolds2bin.tsv ]; then

                echo "${NAME}: running fasta to scaffolds, CONCOCT"
                $dastool/src/Fasta_to_Scaffolds2Bin.sh -i $CONC -e fa \
                        > $DAS/${NAME}_concoct_scaffolds2bin.tsv

        else

                echo "${NAME}: using existing file for CONCOCT: $DAS/${NAME}_concoct_scaffolds2bin.tsv"

        fi

	echo "${NAME}: starting DASTool using ${searcheng}, score threshold ${score}, ${OMP_NUM_THREADS} threads"
        echo "${NAME}: DASTool using $(R --version | head -1)"

	# if DASTool previously failed, remove bin directory
        rm -r $DAS/${NAME}_DASTool_bins

	cd $DAS

        $dastool/DAS_Tool \
                -i $DAS/${NAME}_maxbin2_scaffolds2bin.tsv,$DAS/${NAME}_metabat2_scaffolds2bin.tsv,$DAS/${NAME}_concoct_scaffolds2bin.tsv \
                -l maxbin,metabat,concoct \
                -c $MSPA/contigs.fasta \
                -o ${NAME} \
                --search_engine ${searcheng} \
                --write_bins 1 \
		--write_bin_evals 1 \
                --score_threshold ${score} \
                -t ${OMP_NUM_THREADS}

        # count DASTool bins. If > 0, remove intermediate bins and files
        if [ $(find $DAS/${NAME}_DASTool_bins -name '*.fa' | wc -l) -eq 0 ]; then

                echo "${NAME}: No DASTool bins passing threshold, Aborting"
                qsub -N binning_FAIL_${NAME} -M ${email} ${jobtracker}
		exit 1

        else

                echo "${NAME}: DASTool finished, $(find $DAS/${NAME}_DASTool_bins -name '*.fa' | wc -l) bins passing threshold"

        fi

else

	echo "${NAME}: found $(find $DAS/${NAME}_DASTool_bins -name '*.fa' | wc -l) existing DASTool bins"

fi

############
## CheckM ##
############

echo "${NAME}: checking for existing CheckM2 quality report for DASTool bins"

CHKM=$DAS/CheckM
mkdir -p $CHKM
cd $CHKM

if [ ! -f $CHKM/quality_report.tsv ]; then

	echo "${NAME}: no existing CheckM2 reports found"
	echo "${NAME}: running checkm2 predict"

	# see /workdir/refdbs/checkm2/_README

	source /home/acv46/miniconda3/bin/activate
	conda activate checkm2

	checkm2 predict \
		--input $DAS/${NAME}_DASTool_bins \
		--output-directory $CHKM \
		--allmodels \
		--threads ${OMP_NUM_THREADS} \
		-x fa \
		--force \
		--quiet

	conda deactivate

	echo "${NAME}: checkm2 predict done, see quality_report.tsv"

	if [ ! -f $CHKM/quality_report.tsv ]; then

		echo "${NAME}: ERROR - checkm2 predict failed, Aborting"
                qsub -N CheckM2_predict_FAIL_${NAME} -M ${email} ${jobtracker}
                exit 1

	fi

	echo "${NAME}: checkm2 predict done, see quality_report.tsv"

else

	echo "${NAME}: existing CheckM2 summary files found at ${CHKM}"

fi

echo "${NAME}: checking for existing CheckM abundance profile for DASTool bins"

if [ ! -f $CHKM/abundance_report.tsv ]; then

        echo "${NAME}: no existing CheckM abundance reports found"
	echo "${NAME}: running checkm coverage"

	source /home/acv46/miniconda3/bin/activate
        conda activate checkm

	checkm coverage -q \
                -t ${OMP_NUM_THREADS} \
                -r \
                -x fa \
                -a 0.98 \
                -e 0.02 \
                -m 15 \
                $DAS/${NAME}_DASTool_bins \
                coverage_report.tsv \
                $BAM

	if [ ! -f $CHKM/coverage_report.tsv ]; then

                echo "${NAME}: ERROR - checkm coverage failed, Aborting"
                qsub -N CheckM_coverage_FAIL_${NAME} -M ${email} ${jobtracker}
                exit 1

        fi

	echo "${NAME}: running checkm profile"

	checkm profile -q \
                -f abundance_report.tsv \
                --tab_table \
                coverage_report.tsv

	if [ ! -f $CHKM/abundance_report.tsv ]; then

                echo "${NAME}: ERROR - checkm profile failed, Aborting"
                qsub -N CheckM_profile_FAIL_${NAME} -M ${email} ${jobtracker}
                exit 1

        fi

        echo "${NAME}: checkm profile done, see abundance_report.tsv"

else

	echo "${NAME}: existing CheckM abundance reports found at $CHKM"

fi

#############
## GTDB-Tk ##
#############

# Note: pplacer needs >120 GB memory and at least 4 cores to run

GTD=$DAS/gtdbtk

if [ $gFlag -eq 1 ]; then

	echo "${NAME}: checking for existing GTDB-Tk summary of DASTool bins"

	if [ ! -f $GTD/gtdbtk.bac120.summary.tsv ]; then

		rm -r $GTD
		mkdir -p $GTD
	        cd $GTD

		echo "${NAME}: no existing GTDB-Tk summary found"
		echo "${NAME}: starting GTDB-Tk for DASTool bins"

		export PYTHONPATH=/programs/gtdbtk-2.1.1/lib/python3.6/site-packages/
		export PATH=/programs/gtdbtk-2.1.1/bin:/programs/hmmer/binaries:/programs/prodigal-2.6.3:/programs/FastTree-2.1.10:/programs/fastANI-1.32:/programs/pplacer-Linux-v1.1.alpha19:/programs/mash-Linux64-v2.3:$PATH
		export GTDBTK_DATA_PATH=/workdir/refdbs/gtdbtk/release207_v2/

		# note: the GTDBTK_DATA_PATH directory must be organized to include the following subdirectories:

		##	/workdir/GtdbTK/markers/tigrfam
		##	/workdir/GtdbTK/markers/pfam
		##	/workdir/GtdbTK/msa
		##	/workdir/GtdbTK/masks
		##	/workdir/GtdbTK/pplacer
		##	/workdir/GtdbTK/fastani
		##	/workdir/GtdbTK/taxonomy

		echo "${NAME}: GTDB-Tk - running classify_wf"

		# weird error where multiple cores causes pplacer to crash
		# only when called with qsub / in SGE
		# https://github.com/Ecogenomics/GTDBTk/issues/170#issuecomment-566882748

		gtdbtk classify_wf \
			--cpus ${OMP_NUM_THREADS} \
			--pplacer_cpus ${OMP_NUM_THREADS} \
			--genome_dir $DAS/${NAME}_DASTool_bins \
			--extension fa \
			--out_dir $GTD \

		echo "${NAME}: GTDB-Tk finished, see gtdbtk.bac120.summary.tsv for bin taxa"

	else

		echo "${NAME}: existing GTDB-Tk summary found at $GTD/gtdbtk.bac120.summary.tsv"

	fi

fi

##########
## POEM ##
##########

if [ $pFlag -eq 1 ]; then

	POEM=$DAS/poem
	mkdir -p $POEM
	cd $POEM

	if [ ! -f $POEM/goodbins_${complete}_${contam}.txt ]; then

		echo "${NAME}: subsetting DASTool bins by checkM stats; completeness > ${complete}, contamination < ${contam}"

		grep 'UID' $CHKM/${NAME}_checkM_lineage.txt | \
			awk -v cmp=$complete '$13>cmp' | \
			awk -v cnt=$contam '$14<cnt' | \
			awk '{print $1}' \
			> $POEM/goodbins_${complete}_${contam}.txt

		echo "${NAME}: identified $(wc -l < ${POEM}/goodbins_${complete}_${contam}.txt) bins passing threshold"

		# need to use a diamond version that matches the database build version (0.9.24)
		# reset path to starting path -- otherwise dependency conflicts

		export PATH=$startpath
		export PATH=$HOME/.local/bin:$PATH
		export PYTHONPATH=$HOME/.local/lib/python3.9/site-packages
		source $HOME/miniconda3/bin/activate

		while read line; do

			if [ ! -f $POEM/${line}.fa_output/input.fsa.operon ]; then

				echo "${NAME}: starting POEM operon calling for bin $(grep -n "^${line}$" $POEM/goodbins_${complete}_${contam}.txt | awk -F ":" '{print $1}') of $(wc -l < ${POEM}/goodbins_${complete}_${contam}.txt) -- ${line}.fa"

				#edited script to suppress echo output
				bash /home/acv46/POEM_py3k/bin/run_poem_edit.sh \
					-f $DAS/${NAME}_DASTool_bins/${line}.fa \
					-a n \
					-p pro

				mv $DAS/${NAME}_DASTool_bins/${line}.fa_output $POEM

				echo "${NAME}: POEM operon calling complete for ${line}.fa, see ${line}.fa_output/input.fsa.operon"

			fi

		done < ${POEM}/goodbins_${complete}_${contam}.txt

	fi

fi

##################
## Summary File ##
##################

if [[ $gFlag -eq 1 && ! -f $DAS/${NAME}_BINSTATS.txt ]]; then

	echo "${NAME}: generating summary file of CheckM and GTDB-Tk reports"

	export PATH=/programs/R-4.2.0/bin:$PATH
	Rscript /workdir/users/acv46/stool_PROSeq3/scripts/get_summary.R \
		-q $CHKM/quality_report.tsv \
		-a $CHKM/abundance_report.tsv \
		-g $GTD/gtdbtk.bac120.summary.tsv \
		-n ${NAME} \
		-o $DAS

	if [ ! -f $DAS/${NAME}_BINSTATS.txt ]; then

                echo "${NAME}: ERROR - summary file creation failed, Aborting"
                qsub -N summary_FAIL_${NAME} -M ${email} ${jobtracker}
                exit 1

        fi

        echo "${NAME}: summary file created, $DAS/${NAME}_BINSTATS.txt"


fi

##############
## Clean Up ##
##############

echo "${NAME} -- PIPELINE FINISHED"
echo "${NAME} -- removing intermediate bins and assembly files to save space"
[ ! -d $CONC ] || rm -r $CONC
[ ! -d $BAT ] || rm -r $BAT
[ ! -d $MXB ] || rm -r $MXB

qsub -N CAB_DONE_${NAME} -M ${email} ${jobtracker}


