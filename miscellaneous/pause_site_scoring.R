#!/usr/bin/env Rscript

# see /workdir/users/acv46/EC_PROSeq/pileup/_README
# uses existing merged bams from get_pileup.sh
# subsets bams for MAPQ >= qscore and recounts pileups
# merges pileups into single file

# To do
## fix bind_rows call to get both minus and plus strand results
## incorporate GTF data for genetic context annotation of pause sites
## parallelize apply calls
## rewrite to use single apply call over bin-contig pair instead of two apply calls, for better pbapply reporting

# from Sun 2021 -- "ends of all uniquely mapped RNA reads (bottom lane) were determined and the read count for each 3′ end position was calculated and plotted (top lane). The genomic positions where 3’ end/3’ end median (51-bp window) read counts ratio (pause score) was ≥ 20 and read counts/10^6 reads was ≥ 10 satisfied our stringent definition for a pause site."

# main function to process each contig in the chosen bin set
get_motifs <- function(covdat, bindat, seqdat, gtf, end, range, mincov, minpause, context, method) {
  
  # INPUTS
  ## covdat is coverage data
  ## bindat is all metadata
  ## seqdat is genomic ranges object containing all sequence data
  ## gtf is a gtf annotation file
  ## end is either 3 or 5
  ## range is the context around a position to use in pause score calculation
  ## mincov is the minimum raw coverage value at a position to determine whether a pause score should be calculated
  ## minpause is the the minimum Z score for a position to call it as significantly paused
  ## context is the number of bases surrounding a significant peak
  
  require(tidyverse)
  require(stats)
  require(Biostrings)
  require(pbapply)
  
  # helper function to get pause site z score
  pause_score <- function(row, full, range, con_len, method){
    
    pos <- row[2] %>% as.numeric()
  
    # this code works with both circular and linear sequences
    fullrange <- sort(c(
      c(((pos - 1) - c(1:range)) %% con_len + 1),
      c(((pos - 1) + c(1:range)) %% con_len + 1)
      ))
    
    if (method == "mean") {
    
      mean <- full %>%
        filter(position %in% fullrange) %>%
        select(value) %>% abs() %>% 
        unlist() %>% mean()
    
      std <- full %>%
        filter(position %in% fullrange) %>%
        select(value) %>% abs() %>% 
        unlist() %>% sd()
    
      # calculate z score for pause sites passing raw cov threshold
      zscore <- (abs(as.numeric(row[4])) - mean) / std
    
      list(zscore, as.numeric(row[2]))
    
    } else if (method == "mad") {
  
      mean <- full %>%
        filter(position %in% fullrange) %>%
        select(value) %>% abs() %>% 
        unlist() %>% mean()
    
      mad <- full %>%
        filter(position %in% fullrange) %>%
        select(value) %>% abs() %>% 
        unlist() %>% mad()
    
      # calculate z score for pause sites passing raw cov threshold
      zscore <- (abs(as.numeric(row[4])) - mean) / mad
    
      list(zscore, as.numeric(row[2]))
      
    }
      
  }
  
  # for each bin, apply function to constituent contigs
  pblapply(bindat %>% select(bin) %>% unique() %>% unlist(),
         function(i) {
           pblapply(bindat %>% filter(bin == i) %>% select(contig) %>% unlist(),
                  function(j) {
                    
                    message(paste0("\n","Procesing bin ", i, ", contig ", j))
                    
                    # subset coverage data by contig and sequence end type
                    con_cov <- covdat %>%
                      filter(contig == j & end == end)
                    
                    con_len <- max(con_cov$position) 
                    index <- c((range + 1):(con_len - range - 1))
                    
                    # create separate datasets for plus and minus positions
                    # filter for minimum absolute 
                    # plus_set <- con_cov %>%
                    #   filter(strand == "plus" & abs(value) >= mincov & position %in% index)
                    
                    minus_set <- con_cov %>%
                      filter(strand == "minus" & abs(value) >= mincov & position %in% index)
                    
                    cseq <- seqdat[j][[1]]
                    
                    # get strand-wise z scores for peaks passing raw coverage threshold
                    # eliminate peaks with overlapping ranges by only keeping larger peak
                    ## this is needed to prevent the same sequence from being pulled more than once
                    
                    # plus strand
                    if (nrow(plus_set) > 0) {

                      plus_z <- tibble(z = rep(NA, nrow(plus_set)),
                                       position = rep(NA, nrow(plus_set)))
                      hold <- apply(X = plus_set,
                                    MARGIN = 1,
                                    FUN = pause_score,
                                    full = con_cov %>% filter(strand == "plus"),
                                    range = range,
                                    con_len = con_len,
                                    method = method)

                      plus_z$z <- do.call(rbind, hold)[,1] %>% unlist()
                      plus_z$position <- do.call(rbind, hold)[,2] %>% unlist()
                      plus_z <- plus_z %>%
                        mutate(startrng = position - context) %>%
                        mutate(endrng = position + context) %>%
                        arrange(startrng) %>%
                        filter(z >= minpause)

                      message(paste0("\n","--> plus strand has ", nrow(plus_z), " hits over Z-score threshold"))

                      plus_z$grp <- 1
                      if (nrow(plus_z) > 1) {
                        for (k in 2:nrow(plus_z)) {

                          if (plus_z$endrng[k - 1] >= plus_z$startrng[k]) {

                            plus_z$grp[k] <- plus_z$grp[k - 1]

                          } else {

                            plus_z$grp[k] <- plus_z$grp[k - 1] + 1

                          }
                        }
                      }

                      message(paste0("\n","----> plus strand has ", nrow(plus_z), " hits remaining after overlap removal"))

                      if (nrow(plus_z) > 0) {
                        plus_z <- plus_z %>%
                          group_by(grp) %>%
                          top_n(1, z) %>%
                          ungroup() %>%
                          mutate(bin = i) %>%
                          mutate(contig = j) %>%
                          mutate(strand = "plus")
                          plus_z$sequence <- apply(X = plus_z, MARGIN = 1,
                                                   FUN = function(x) {
                                                     cseq[x[3]:x[4]] %>%
                                                       as.character()
                                                     }
                                                   )

                        plus_z

                      }

                    } %>% bind_rows()
                    
                    # minus strand
                    if (nrow(minus_set) > 0) {

                      minus_z <- tibble(z = rep(NA, nrow(minus_set)),
                                        position = rep(NA, nrow(minus_set)))
                      hold <- apply(X = minus_set,
                                    MARGIN = 1,
                                    FUN = pause_score,
                                    full = con_cov %>% filter(strand == "minus"),
                                    range = range,
                                    con_len = con_len,
                                    method = method)

                      minus_z$z <- do.call(rbind, hold)[,1] %>% unlist()
                      minus_z$position <- do.call(rbind, hold)[,2] %>% unlist()
                      minus_z <- minus_z %>%
                        mutate(startrng = position - context) %>%
                        mutate(endrng = position + context) %>%
                        arrange(startrng) %>%
                        filter(z >= minpause)

                      message(paste0("\n","--> minus strand has ", nrow(minus_z), " hits over Z-score threshold"))

                      minus_z$grp <- 1
                      if (nrow(minus_z) > 1) {
                        for (k in 2:nrow(minus_z)) {

                          if (minus_z$endrng[k - 1] >= minus_z$startrng[k]) {

                            minus_z$grp[k] <- minus_z$grp[k - 1]

                          } else {

                            minus_z$grp[k] <- minus_z$grp[k - 1] + 1

                          }
                        }
                      }

                      message(paste0("\n","----> minus strand has ", nrow(minus_z), " hits remaining after overlap removal"))

                      if (nrow(minus_z) > 0) {
                        minus_z <- minus_z %>%
                          group_by(grp) %>%
                          top_n(1, z) %>%
                          ungroup() %>%
                          mutate(bin = i) %>%
                          mutate(contig = j) %>%
                          mutate(strand = "minus")
                          minus_z$sequence <- apply(X = minus_z, MARGIN = 1,
                                                   FUN = function(x) {
                                                     cseq[x[3]:x[4]] %>%
                                                       reverseComplement() %>%
                                                       as.character()
                                                     }
                                                   )

                        minus_z

                      }

                    } %>% bind_rows()
                    
                }) %>% bind_rows()
       
        }) %>% bind_rows()
  
  
}

