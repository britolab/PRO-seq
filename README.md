## PRO-seq and RNAseq in *E. coli* and microbiomes
Much of the code in this repo was not used in the final analyses. The most important scripts are listed below.

### metagenome assembly

[`run_clean_assemble_bin.sh`](https://github.com/britolab/PRO-seq/blob/main/metagenome_assembly/run_clean_assemble_bin.sh) is the master script for read QC, metagenome assembly, contig binning, bin QC, and taxonomic assignment. Parts of this script are hard-coded to work with the Cornell BioHPC SGE scheduler and the Brito Lab server structure.

### bin annotation

Genes were called in metagenomic bins using [`run_prokka.sh`](https://github.com/britolab/PRO-seq/blob/main/miscellaneous/run_prokka.sh). gtf annotations output by prokka can be converted to R objects using [`gtf2tibble.R`](https://gist.github.com/acvill/03343034392cff158d2369483ed8935f).

### transcript mapping

The main scripts for cleaning up transcript reads and mapping those reads to references can be found in the [Danko Lab proseq2.0 repo](https://github.com/Danko-Lab/proseq2.0). Once you have bam files, per-base coverage reports can be generated with [`get_pileup_correct.sh`](https://github.com/britolab/PRO-seq/blob/main/miscellaneous/get_pileup_correct.sh).

### data processing and visualization

[`EC_peaks.rmd`](https://github.com/britolab/PRO-seq/blob/main/data_processing_and_figures/EC_peaks.rmd) and [`Stool_PRO-seq.Rmd`](https://github.com/britolab/PRO-seq/blob/main/data_processing_and_figures/Stool_PRO-seq.Rmd) contain the R code used for the *E. coli* and human microbiome analyses, respectively. The Rmarkdown documents are ordered by main sections (`#`) and subsections (`##`/`###`).

### data availability

*E. coli* sequencing reads: https://www.ncbi.nlm.nih.gov/sra/PRJNA800038  
microbiome sequencing reads: https://www.ncbi.nlm.nih.gov/sra/PRJNA800070
