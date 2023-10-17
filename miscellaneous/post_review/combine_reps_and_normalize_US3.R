setwd("/workdir/users/acv46/stool_PROSeq3/transcriptomes")
library(tidyverse)

read_cov <- function(file_) {

  base_ <- read_tsv(file = file_, col_names = TRUE)
  PRO_sum  <- sum(base_[,grepl("PRO", names(base_))])
  TEXp_sum <- sum(base_[,grepl("TEXp", names(base_))])
  TEXm_sum <- sum(base_[,grepl("TEXm", names(base_))])
  
  base_ |>
    mutate(PRO_f = (PRO_A_f + PRO_B_f + PRO_C_f) / PRO_sum,
           PRO_r = (PRO_A_r + PRO_B_r + PRO_C_r) / PRO_sum,
           TEXp_f = (TEXp_A_f + TEXp_B_f + TEXp_C_f) / TEXp_sum,
           TEXp_r = (TEXp_A_r + TEXp_B_r + TEXp_C_r) / TEXp_sum,
           TEXm_f = (TEXm_A_f + TEXm_B_f + TEXm_C_f) / TEXm_sum,
           TEXm_r = (TEXm_A_r + TEXm_B_r + TEXm_C_r) / TEXm_sum) |>
    select(contig, position, PRO_f, PRO_r,
           TEXp_f, TEXp_r, TEXm_f, TEXm_r) |>
    pivot_longer(cols = -c("contig","position"), names_to = "lib", values_to = "value") |>
    separate(col = lib, into = c("type","strand"), sep = "_") |>
    mutate(strand = ifelse(strand == "f", "+", "-"))
  
}

write_tsv(x = read_cov(file_ = "US3/coverage/US3_combined_depth.txt"),
          file = "US3/coverage/US3_combined_depth_normalized.txt",
          col_names = TRUE)