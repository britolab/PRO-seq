#$ -S /bin/bash
#$ -N graphlan
#$ -V
#$ -e /workdir/users/acv46/stool_PROSeq2/log/graphlan_$JOB_ID.err
#$ -o /workdir/users/acv46/stool_PROSeq2/log/graphlan_$JOB_ID.out
#$ -wd /workdir/users/acv46/stool_PROSeq2
#$ -l h_vmem=80G
#$ -q long.q@cbsubrito
#$ -t 1

# input metaphlan2 profiles (run_metaphlan.sh)
# output hclust2 heatmaps and graphlan cladogram

#for i in /workdir/users/acv46/stool_PROSeq2/out/*rep*/metaphlan/*_profile.txt; do
	#new=$(echo $i | sed 's|_mgm||g')
	#mv $i $new
#done

MPH=/workdir/users/acv46/stool_PROSeq2/out/graphlan
mkdir -p $MPH

# hard-coded, copy metaphlan profiles from existing directories
############################
cp /workdir/users/acv46/stool_PROSeq2/out/US*/metaphlan/*_profile.txt $MPH
cp /workdir/users/acv46/mgmAssembly/metaphlan/US3_3Nov2020_mgm_profile.txt $MPH
cp /workdir/users/acv46/mgmAssembly/metaphlan/US2_5Nov2020_mgm_profile.txt $MPH
############################

echo "checking for metaphlan profiles at ${MPH}"
echo "$(ls $MPH | grep "profile.txt" | wc -l) metaphlan profiles found"
echo "merging metaphlan profiles"

export PATH=/programs/MetaPhlAn-2.0:/programs/MetaPhlAn-2.0/utils:$PATH
export OMP_NUM_THREADS=8

merge_metaphlan_tables.py ${MPH}/*_profile.txt | \
	sed "s/_profile//g" > ${MPH}/merged_abundance_table.txt

echo "formatting abundance tables for hclust input"

# species-level
grep -E "(s__)|(^ID)" ${MPH}/merged_abundance_table.txt | \
	grep -v "t__" | \
	sed 's/^.*s__//g' > ${MPH}/merged_abundance_table_species.txt

# genus-level
grep -E "(g__)|(^ID)" ${MPH}/merged_abundance_table.txt | \
	grep -v "s__" | \
	sed 's/^.*g__//g' > ${MPH}/merged_abundance_table_genus.txt

# family-level
grep -E "(f__)|(^ID)" ${MPH}/merged_abundance_table.txt | \
        grep -v "g__" | \
        sed 's/^.*f__//g' > ${MPH}/merged_abundance_table_family.txt

# phylum-level
grep -E "(p__)|(^ID)" ${MPH}/merged_abundance_table.txt | \
        grep -v "c__" | \
        sed 's/^.*p__//g' > ${MPH}/merged_abundance_table_phylum.txt

echo "creating hclust heatmaps"

source $HOME/miniconda3/bin/activate

# if needed, make local installations
#conda install -c biobakery hclust2

for table in ${MPH}/merged_abundance_table_*.txt; do

	taxa=$(basename $table | sed "s/merged_abundance_table_//g" | sed "s/.txt//g")

	hclust2.py \
		-i ${table} \
		-o ${MPH}/abundance_heatmap_${taxa}.png \
		--ftop 25 \
		--f_dist_f braycurtis \
		--s_dist_f braycurtis \
		--cell_aspect_ratio 0.5 \
		-l \
		--flabel_size 6 \
		--slabel_size 6 \
		--max_flabel_len 100 \
		--max_slabel_len 100 \
		--minv 0.1 \
		--dpi 300

	echo "hclust heatmap created for taxonomic level = ${taxa}"

done

echo "creating graphlan input files"
# note: conda install of graphlan won't resolve
#pip install --user export2graphlan
#pip install --user hclust2

E2G=/home/britolab/acv46/.local/bin/export2graphlan.py
source /programs/miniconda3/bin/activate graphlan

# annotations use taxonomic levels
# 2 = phylum, 7 = species
# metadata rows are 0-indexed

python2.7 $E2G \
	-i ${MPH}/merged_abundance_table.txt \
	--tree ${MPH}/merged_abundance.tree.txt \
	--annotation ${MPH}/merged_abundance.annot.txt \
	--most_abundant 100 \
	--abundance_threshold 1 \
	--least_biomarkers 10 \
	--annotations 6 \
	--external_annotations 7 \
	--background_levels 3 \
	--min_clade_size 1 \
	--metadata_rows 0 \
	--skip_rows 1

echo "annotating graphlan cladogram"

graphlan_annotate.py \
	--annot ${MPH}/merged_abundance.annot.txt \
	${MPH}/merged_abundance.tree.txt \
	${MPH}/merged_abundance.xml

echo "printing graphlan cladogram"

graphlan.py \
	--dpi 300 \
	${MPH}/merged_abundance.xml \
	${MPH}/merged_abundance.png \
	--external_legends
