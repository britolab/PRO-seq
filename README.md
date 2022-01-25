## PRO-seq and RNAseq in *E. coli* and microbiomes
Much of the code in this repo was not used in the final analyses. The most important scripts are list below.

### metagenome assembly

[`run_clean_assemble_bin.sh`](https://github.com/britolab/PRO-seq/blob/main/metagenome_assembly/run_clean_assemble_bin.sh) is the master script for read QC, metagenome assembly, contig binning, and bin QC. Parts of this script are hard-coded to work with the Cornell BioHPC SGE scheduler and the Brito Lab server structure.

### data 

*E. coli* sequencing reads: https://www.ncbi.nlm.nih.gov/sra/PRJNA800038  
microbiome sequencing reads: https://www.ncbi.nlm.nih.gov/sra/PRJNA800070
