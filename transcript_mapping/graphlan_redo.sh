#!/bin/bash

# to be used in conjunction with files made in run_graphlan.sh
# usage: sh graphlan_redo.sh suffix

MPH=/workdir/users/acv46/stool_PROSeq2/out/graphlan
mkdir -p $MPH/temp

E2G=/home/britolab/acv46/.local/bin/export2graphlan.py
source /programs/miniconda3/bin/activate graphlan

# annotations use taxonomic levels
# 2 = phylum, 7 = species
# metadata rows are 0-indexed

python2.7 $E2G \
        -i ${MPH}/merged_abundance_table.txt \
        --tree ${MPH}/temp/merged_abundance.tree_${1}.txt \
        --annotation ${MPH}/temp/merged_abundance.annot_${1}.txt \
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
        --annot ${MPH}/temp/merged_abundance.annot_${1}.txt \
        ${MPH}/temp/merged_abundance.tree_${1}.txt \
        ${MPH}/temp/merged_abundance_${1}.xml

echo "printing graphlan cladogram"

graphlan.py \
        --dpi 300 \
        ${MPH}/temp/merged_abundance_${1}.xml \
        ${MPH}/figs/merged_abundance_${1}.png \
        --external_legends

#rm -r ${MPH}/temp
