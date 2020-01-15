#!/usr/bin/env Rscript

Args <- commandArgs(TRUE)
rounds <- as.numeric(Args[1])

source("./functions/functions_reports.R")
options(bitmapType='cairo')

r_ <- (rounds %% 10) + 1

ed_ <- (rounds - 1) %/% 10

pat.tmp <- NULL

if(ed_ == 0) pat.tmp <- "_g-0.3k-5.2k"
if(ed_ == 1) pat.tmp <- "_g-0.3k-4.68k"
if(ed_ == 2) pat.tmp <- "_g-0.3k-5.72k"
if(ed_ == 3) pat.tmp <- "_g-0.3k-6.24k"

pattern <- paste0("r-", r_, pat.tmp)

sampled.path <- "./data/"
sampled.names <- list.files(path = sampled.path, pattern = "*.RData")
r.this <- sampled.names[grepl(pattern, sampled.names)]
r.this <- r.this[grepl("_g-", r.this)]
r.this <- r.this[!grepl("coefs.wb.", r.this)]

t <- Sys.time()
for(sampled in r.this){
  load(paste(sampled.path, sampled, sep = ""))
  brain_case@name %>% paste0("'s network is being plotted") %>% print()
  reports_netviz(brain_case)
}

paste("Netviz of round", r_,
      "of", pat.tmp,
      "index", rounds,
      "took", (Sys.time()-t)) %>% print()
