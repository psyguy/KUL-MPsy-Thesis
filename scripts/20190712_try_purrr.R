rm(list = ls())

# setwd("I:/Thesis_Codes/Thesis_R_Codes/")

source("./functions/functions_extract.R")


path.to.brains <- "./data/5200-edges"
path.to.snp <- "./data/5200-edges-snapshots"


pattern <- "_g-0.3k-5.2k"

brain_files <- list.files(path = path.to.brains,
                          pattern = "*.RData")

r.this <- brain_files[grepl(pattern, brain_files)]
r.this <- r.this[grepl("_g-", r.this)]
# r.this <- r.this[grepl(paste0("_r-", r_), r.this)]
r.this <- r.this[!grepl("coefs.wb.", r.this)]
brain_locations <- path.to.brains %>% 
  paste(r.this, sep = "/")

# Sys.time() %>% print()

for(i in 1:25){
  ind <- ((i-1)*20+1):(i*20)
  paste(ind[[1]], "to", ind[[20]], "is being read now") %>% print()
  this.brain_location <- brain_locations[ind]
  
  Sys.time() %>% print()
  system.time(snp <- this.brain_location %>%
                ldply(extract_brains,
                      snapshots =  seq(5e3, 100e3, 5e3)
                )
              ) %>% print()
  
  save_vars(
    list.of.vars = "snp",
    prefix = paste0("snp_", ind[[1]], "-", ind[[20]]),
    path = path.to.snp
  )
}
