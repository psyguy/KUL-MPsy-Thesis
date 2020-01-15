#!/usr/bin/env Rscript

# Args <- commandArgs(TRUE)
# rounds <- as.numeric(Args[1])
rm(list = ls())

source("./functions/functions_reports.R")

pattern <- "coefs.wb.hpc_"

paste("Now reading round", r_, "of", pat.tmp, "index", rounds) %>% print()

 sampled.path = "./data/harvest_hpc/"

sampled.names <- list.files(path = sampled.path, pattern = "*.RData")
r.this <- sampled.names[grepl(pattern, sampled.names)]
  # r.this <- r.this[grepl("_g-", r.this)]
  # r.this <- r.this[!grepl("coefs.wb.", r.this)]
  
coefs.wb.all <- NULL
t <- Sys.time()
for(sampled in r.this){
  load(paste(sampled.path, sampled, sep = ""))
  coefs.wb.all <- coefs.wb.all %>% rbind(coefs.wb)
}
(Sys.time()-t) %>% paste("for", pattern) %>% print()


# making the verbal description of the parameter sets
coefs.all <- coefs.wb.all %>% mutate(`Verbal Description` = 
                                       paste(`Epsilon Proportion`, `a Proportion`, sep="_"))
coefs.all$`Verbal Description` <- coefs.all$`Verbal Description` %>% as.factor()

levels(coefs.all$`Verbal Description`) <- c("Hyper-coupled Minority",
                                            "Hyper-chaotic Minority",
                                            "Homogeneous Society",
                                            "Hypo-chaotic Minority",
                                            "Hypo-coupled Minority")

save_vars("coefs.all", prefix = "coefs.all")
