rm(list = ls())
load("./data-pc/5200-edges/life-10_eps-1v5v0_a-0v6v0_r-9_g-0.3k-5.2k_Hannah Bauwens_20190709_2225.RData")

# setwd("I:/Thesis_Codes/Thesis_R_Codes/")

source("./functions/functions_extract.R")


path.to.brains <- "./data-pc/5200-edges"
path.to.snp <- "./data/5200-edges-snapshots"
pattern <- "_g-0.3k-5.2k"
brain_files <- list.files(path = path.to.brains,
                          pattern = "*.RData")
brain_files_1e6 <- brain_files[grepl("life-10", brain_files)]

r.this <- brain_files[grepl("life-10", brain_files)]
r.this <- r.this[!grepl("coefs.wb.", r.this)]
brain_locations <- path.to.brains %>%
  paste(r.this, sep = "/")
path.to.snp <- "./data-pc/1e6-snapshots"
brain_files <- list.files(path = path.to.brains,+pattern = "*.RData")
r.this <- brain_files[grepl("life-10", brain_files)]
brain_locations <- path.to.brains %>%
  paste(r.this, sep = "/")
b <- brain_case




extract_now <- function(b){
  iden <- extract_identifiers(b)$whole
  owner <- iden$Owner
  verbal <- iden$Verbal.Description
  
  m <- b@now$mat.connectivity
  a <- b@now$activities
  
  o <- list(
            owner = owner,
            verbal = verbal,
            m = m,
            a = a
            )
  o %>% return()
}

dist_maker <- function(l.now){
  num.persons <- 50
  dist.m <- matrix(0, num.persons, nrow = num.persons)
  dist.a <- matrix(0, num.persons, nrow = num.persons)
  system.time(
    for(i in 1:num.persons){
      m1 <- s[[i]] %>% my_seriate()
      print(i)
      for(j in i:(num.persons)){
        m2 <- s[[j]] %>% my_seriate()
        dist.m[j,i] <- 1 - MatrixCorrelation::PSI(m1,m2)
        dist.a[j,i] <- 1 - MatrixCorrelation::PSI(m1,m2)
        print(j)
      }
    }
  )
  
}

o <- extract_now(b)

for (i in 1:25) {
  ind <- ((i - 1) * 20 + 1):(i * 20)
  paste(ind[[1]], "to", ind[[20]], "is being read now") %>% print()
  this.brain_location <- brain_locations[ind]
  
  Sys.time() %>% print()
  system.time(snp <- this.brain_location %>%
                ldply(extract_brains,
                      snapshots =  seq(5e3, 100e3, 5e3))) %>% print()
  
  save_vars(
    list.of.vars = "snp",
    prefix = paste0("snp_", ind[[1]], "-", ind[[20]]),
    path = path.to.snp
  )
}
