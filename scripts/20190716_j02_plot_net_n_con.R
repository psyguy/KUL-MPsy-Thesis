#!/usr/bin/env Rscript
# 
# Args <- commandArgs(TRUE)
# ind <- as.numeric(Args[1])

rm(list = ls())

source("./functions/functions_extract.R")

path.to.snapshots <- "./data/5200-edges-snapshots/"

snp_files <- list.files(path = path.to.snapshots, pattern = "*.RData")


# snp_files <- snp_files[ind]

t0 <- Sys.time()
for(this.snp in snp_files[6:10]){
  this.snp %>% print()
  load(paste0(path.to.snapshots, this.snp))

  snp.whole <- snp %>%
    filter(Partition == "whole")
  rew.un <- snp.whole$Rewiring %>% unique() %>% sort()
  
  t <- Sys.time()
  for(i in 1:nrow(snp.whole)){
    if(snp.whole$Rewiring[[i]] %% 25e3) next
    owner <- snp.whole$Owner[[i]] %>% as.character()
    vd <- snp.whole$Verbal.Description[[i]] %>% as.character()
    rew <- snp.whole$Rewiring[[i]] / 1000
    rew <- rew %>% str_pad(width = 4, pad = "0")
    v <- snp.whole$adj.mat.vect[[i]]
    
    title <- paste0(owner,
                    " (", tolower(vd), ")",
                    " at ", rew, "k rewirings")
    extract_plotcon(v, title = title, path.fig = "figures/connectivity-plots_25k")
    extract_plotnet(v, title = title, path.fig = "figures/connectivity-plots_25k")
  }
  (Sys.time()-t) %>% print()
  
  
}
(Sys.time()-t0) %>% print()
