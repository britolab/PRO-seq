#!/usr/bin/Rscript

library(tidyverse)

setwd("/workdir/users/acv46/stool_PROSeq2/deseq/bam/pileup")

tRNA_all_goodbins <- read_csv(file = "tRNA_all_goodbins_534.txt", col_names = T)

message(paste0("total tRNA file read in"))

US2_contigs_tRNA <- tRNA_all_goodbins %>%
	filter(samp == "US2") %>%
	select(contig) %>% unique %>% as_vector
US3_contigs_tRNA <- tRNA_all_goodbins %>%
	filter(samp == "US3") %>%
	select(contig) %>% unique %>% as_vector

message(paste0("individual tRNA objects created"))

read_csv(file = "US2_coverage_q30_correct2_processed.txt", 
		col_names = T) %>% 
	filter(name %in% samps) %>%
	filter(contig %in% US2_contigs_tRNA) %>%
	mutate(nlab = ifelse(end == "full",
			paste0(type),
			paste0(type, " ", end, "' end"))) %>%
	write_csv(file = "US2_coverage_q30_correct2_processed_tRNAs.txt",
			col_names = T)

message(paste0("US2 subset written to US2_coverage_q30_correct2_processed_tRNAs.txt"))

read_csv(file = "US3_coverage_q30_correct2_processed.txt",
                col_names = T) %>%
        filter(name %in% samps) %>%
        filter(contig %in% US3_contigs_tRNA) %>%
        mutate(nlab = ifelse(end == "full",
                        paste0(type),
                        paste0(type, " ", end, "' end"))) %>%
        write_csv(file = "US3_coverage_q30_correct2_processed_tRNAs.txt",
                        col_names = T)

message(paste0("US3 subset written to US3_coverage_q30_correct2_processed_tRNAs.txt"))
