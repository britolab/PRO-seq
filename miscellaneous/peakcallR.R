#!/usr/bin/env Rscript

# To do
## fix bind_rows call to get both minus and plus strand results
## incorporate GTF data for genetic context annotation of pause sites
## parallelize apply calls
## rewrite to use single apply call over bin-contig pair instead of two apply calls, for better pbapply reporting

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(stats))
suppressMessages(library(Biostrings))
#library(pbapply)

option_list=list(
  make_option(c("-o", "--out"),
              action = "store",
              type = "character",
              default = NULL,
              help = "REQUIRED: name for csv results file ({out}.csv)"),
  make_option(c("-d", "--outdir"),
              action = "store",
              type = "character",
              default = NULL,
              help = "REQUIRED: path to write results. Will create the directory if it does not already exist."),
  make_option(c("-c", "--covdat"),
              action = "store",
              type = "character",
              default = NULL,
              help = "REQUIRED: per-base coverage file (csv). Give full path if file not in outdir."),
  make_option(c("-b", "--bindat"),
              action = "store",
              type = "character",
              default = NULL,
              help = "REQUIRED: bin-to-contig mapping file and metadata (csv). Give full path if file not in outdir."),
  make_option(c("-f", "--fasta"),
              action = "store",
              type = "character",
              default = NULL,
              help = "REQUIRED: fasta, deflines must match bin data contigs. Give full path if file not in outdir."),
  make_option(c("-e", "--end"),
              action = "store",
              type = "integer",
              default = 3,
              help = "either 3 or 5, representing which read end to use for peak calls (default: 3)"),
  make_option(c("-r", "--range"),
              action = "store",
              type = "integer",
              default = 25,
              help = "the sequence context surrounding a peak used to calculate the score (default: 25)"),
  make_option(c("-y", "--mincov"),
              action = "store",
              type = "integer",
              default = 10,
              help = "the minimum normalized coverage at a position for it to be considered for peak calling (default: 10)"),
  make_option(c("-z", "--minscore"),
              action = "store",
              type = "integer",
              default = 5,
              help = "the minimum z-score for a peak to be included in the output (default: 5)"),
  make_option(c("-s", "--seqcontext"),
              action = "store",
              type = "integer",
              default = 20,
              help = "the nt sequence context around each peak to include in the output (default: 20)"),
  make_option(c("-m", "--method"),
              action = "store",
              type = "character",
              default = "modz",
              help = "the method to use for peak score calculation, either 'z' for z-score or 'modz' MAD-modified z-score (default: 'modz')"))
opt_parser=OptionParser(option_list = option_list)
opt=parse_args(opt_parser)

# main function to process each contig in the chosen bin set
get_motifs <- function(covdat, bindat, fasta, end, range,
                       mincov, minscore, context, method) {
  
  # helper function to get pause site z score
  pause_score <- function(row, full, range, con_len, method){
    
    pos <- row[2] %>% as.numeric()
  
    # this code works with both circular and linear 1-based sequence indices
    fullrange <- sort(c(
      c(((pos - 1) - c(1:range)) %% con_len + 1),
      c(((pos - 1) + c(1:range)) %% con_len + 1)
      ))
    
    flank <- full %>%
      filter(position %in% fullrange) %>%
      select(value) %>% abs() %>% unlist()
    
    if (method == "z") {
    
      Mn <- base::mean(flank)
      St <- stats::sd(flank)
      zscore <- (abs(as.numeric(row[4])) - Mn) / St
    
      list(round(zscore,3), as.numeric(row[2]))
    
    } else if (method == "modz") {
  
      # note: mad() calculates 1.4826 * median(abs(dat - med))
      # this method will give z = Inf if more than half of values in range are equal
      Md <- stats::median(flank)
      MAD <- stats::mad(flank)
      zscore <- (abs(as.numeric(row[4])) - Md) / MAD
    
      list(round(zscore,3), as.numeric(row[2]))
      
    }
      
  }
  
  
  bindat <- read_csv(file = bindat, col_names = T, show_col_types = F)
  covdat <- read_csv(file = covdat, col_names = T, show_col_types = F)
  seqdat <- readDNAStringSet(filepath = fasta, format = "fasta")
  #totres <- tibble()
  
  # for each bin, apply function to constituent contigs
  as_tibble(do.call(bind_rows, 
  lapply(bindat %>% select(contig) %>% unique() %>% unlist(),
           function(j) {
           
             i <- bindat %>% filter(contig == j) %>%
               select(bin) %>% unlist()
             
             message(paste0("Processing bin = ", i, ", contig = ", j))
                      
             # subset coverage data by contig and sequence end type
             con_cov <- covdat %>%
               filter(contig == j & end == end)
             con_len <- max(con_cov$position) 
             index <- c((range + 1):(con_len - range - 1))
            
             # create separate datasets for plus and minus positions
             # filter for minimum absolute coverage
             plus_set <- con_cov %>%
               filter(strand == "plus" & abs(value) >= mincov & position %in% index)
            
             minus_set <- con_cov %>%
               filter(strand == "minus" & abs(value) >= mincov & position %in% index)
             
             cseq <- seqdat[j][[1]]
             
             results <- tibble()
             
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
                 filter(z >= minscore)
  
               message(paste0("--> plus strand has ", nrow(plus_z), " hits over Z-score threshold"))
  
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
  
               message(paste0("----> plus strand has ", nrow(plus_z), " hits remaining after overlap removal"))
  
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
                   
                   results <- bind_rows(plus_z, results)
                   
               }
  
             }
            
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
                 filter(z >= minscore)
  
               message(paste0("--> minus strand has ", nrow(minus_z), " hits over Z-score threshold"))
  
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
  
               message(paste0("----> minus strand has ", nrow(minus_z), " hits remaining after overlap removal"))
  
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
                  
                  results <- bind_rows(minus_z, results)
  
               } 
  
             }
             
             results
             
           })))
  
}

writeout <- gsub("$/","",opt$outdir)
if (!dir.exists(writeout)) {
  dir.create(writeout)
}

setwd(writeout)

# run function on input data
get_motifs(covdat   = opt$covdat,
           bindat   = opt$bindat,
           fasta    = opt$fasta,
           end      = opt$end,
           range    = opt$range,
           mincov   = opt$mincov,
           minscore = opt$minscore,
           context  = opt$seqcontext,
           method   = opt$method) %>% 
  write_csv(., file = paste0(opt$out,".csv"), col_names = T)

message(paste0("\n","Done! Results written to ", writeout,"/",opt$out,".csv"))