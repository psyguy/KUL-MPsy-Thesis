#!/usr/bin/env Rscript

# Args <- commandArgs(TRUE)
# rounds <- as.numeric(Args[1])



source("./functions/functions_reports.R")


r_ <- (rounds %% 10) + 1

ed_ <- (rounds - 1) %/% 10

pat.tmp <- NULL

if(ed_ == 0) pat.tmp <- "_g-0.3k-5.2k"
if(ed_ == 1) pat.tmp <- "_g-0.3k-4.68k"
if(ed_ == 2) pat.tmp <- "_g-0.3k-5.72k"
if(ed_ == 3) pat.tmp <- "_g-0.3k-6.24k"

pattern <- paste0("r-", r_, pat.tmp)

paste("Now reading round", r_, "of", pat.tmp, "index", rounds) %>% print()

 sampled.path = "./data/"

  sampled.names <- list.files(path = sampled.path, pattern = "*.RData")
  r.this <- sampled.names[grepl(pattern, sampled.names)]
  r.this <- r.this[grepl("_g-", r.this)]
  r.this <- r.this[!grepl("coefs.wb.", r.this)]
  
  coefs.wb <- NULL
  t <- Sys.time()
  for(sampled in r.this){
    # rm(brain_case)
    load(paste(sampled.path, sampled, sep = ""))
    brain_case@name %>% print()
      for(i in 1:(brain_case@age$rewires-1)){
        m <- brain_case@history$mat.connectivity[[i]]
        if(is.null(m)) next
        if(!(i%%5000)) paste("Adding coefs of rewiring", i,
                             "for", brain_case@name) %>% print()
        coefs.wb <- coefs.wb %>% rbind(m %>% netmeas_wbcoefs(parameters = brain_case@parameters,
                                                                 name = brain_case@name,
                                                                 rewiring = i))
      }
  }
  (Sys.time()-t) %>% paste("for", pattern) %>% print()

  save_vars("coefs.wb", prefix = paste0("coefs.wb.hpc_", pattern))
