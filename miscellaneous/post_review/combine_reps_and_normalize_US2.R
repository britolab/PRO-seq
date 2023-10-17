library(tidyverse)

read_cov <- function(file_, outdir_, nchunks_) {

  fname_ <- gsub(".txt", "", basename(file_))
  base_ <- read_tsv(file = file_, col_names = TRUE)
  PRO_sum  <- sum(base_[,grepl("PRO", names(base_))])
  TEXp_sum <- sum(base_[,grepl("TEXp", names(base_))])
  TEXm_sum <- sum(base_[,grepl("TEXm", names(base_))])
  rows_ <- nrow(base_)
  
  # remove original and read in again
  # too much memory
  rm(base_)
  
  # read in again, split into n chunks and process separately
  
  split_ <- split(x = read_tsv(file = file_, col_names = TRUE),
                  f = ceiling((1:rows_ / rows_) * nchunks_))
  
  for (i in 1:nchunks_) {
    
    message(paste("Processing chunk", i, "of", nchunks_))
    chunk_ <- sprintf(paste0("%0",nchar(nchunks_),"d"),i)
    
    mutate(split_[[i]],
           PRO_f = (PRO_A_f + PRO_B_f + PRO_C_f) / PRO_sum * 1e9,
           PRO_r = (PRO_A_r + PRO_B_r + PRO_C_r) / PRO_sum * 1e9,
           TEXp_f = (TEXp_A_f + TEXp_B_f + TEXp_C_f) / TEXp_sum * 1e9,
           TEXp_r = (TEXp_A_r + TEXp_B_r + TEXp_C_r) / TEXp_sum * 1e9,
           TEXm_f = (TEXm_A_f + TEXm_B_f + TEXm_C_f) / TEXm_sum * 1e9,
           TEXm_r = (TEXm_A_r + TEXm_B_r + TEXm_C_r) / TEXm_sum * 1e9) |>
      select(contig, position, PRO_f, PRO_r,
             TEXp_f, TEXp_r, TEXm_f, TEXm_r) |>
      pivot_longer(cols = -c("contig","position"), names_to = "lib", values_to = "value") |>
      separate(col = lib, into = c("type","strand"), sep = "_") |>
      mutate(strand = ifelse(strand == "f", "+", "-")) |>
      write_tsv(file = paste0(outdir_, fname_, "_chunk", chunk_, ".txt"))
    
  }
  
}

setwd("/workdir/users/acv46/stool_PROSeq3/transcriptomes")
read_cov(file_ = "US2/coverage/US2_combined_depth.txt", outdir_ = "US2/coverage/", nchunks_ = 10)
read_cov(file_ = "US3/coverage/US3_combined_depth.txt", outdir_ = "US3/coverage/", nchunks_ = 10)