#!/usr/bin/env Rscript

Args <- commandArgs(TRUE)
counter <- as.numeric(as.character(Args[1]))


# rm(list=ls())

library(tidyverse)
library(dplyr)
library(plyr)
source("./functions/functions_save_vars.R")

library(HHG)

# scaling feature vectors -------------------------------------------------

## do not run this chunck on HPC, already calculated

load("./data/StatusQuo_20191129_0425.RData")
scale.features <- function(features, dont.scale = TRUE){
  if(dont.scale) features %>% return()
  f <- do.call(rbind, features)
  f.center <- f %>% apply(2, mean)
  f.sd <- f %>% apply(2, sd)

  features.scaled <- features %>% map(scale, f.center, f.sd)
  features.scaled %>% return()
}

features.connectivities <- features.connectivities %>% scale.features()
features.activities <- features.activities %>% scale.features()

distances.connectivities <- features.connectivities %>%
  map(~as.matrix(dist(.x, diag = TRUE, upper = TRUE)))
distances.activities <- features.activities %>%
  map(~as.matrix(dist(.x, diag = TRUE, upper = TRUE)))

d.c.unscaled <- distances.connectivities
d.a.unscaled <- distances.activities

save_vars(c("d.c.unscaled", "d.a.unscaled"),
          prefix = "distances of unscaled features"
          )


# calculating HHG per pair ------------------------------------------

load("./data/distances of unscaled features_20191215_1419.RData")

i.and.j <- list()
k <- 1
for(i in 1:length(d.c)){
  for(j in (i):length(d.c)){
    i.and.j[[k]] <- c(i,j)
    k <- k + 1
  }
}

i <- i.and.j[[counter]][1]
j <- i.and.j[[counter]][2]

system.time(
hhg.con <- list(indices = c(i,j),
                names = c(names(d.c[i]),names(d.c[j])),
                hhg.values = hhg.test(d.c[[i]], d.c[[j]], nr.perm = 20000)
                )
)
system.time(
hhg.act <- list(indices = c(i,j),
                names = c(names(d.a[i]),names(d.a[j])),
                hhg.values = hhg.test(d.a[[i]], d.a[[j]], nr.perm = 20000)
                )
)
save_vars(c("hhg.con", "hhg.act"),
          prefix = paste0("HHG-unscaled_", counter, "_indices-", i, "-", j)
          )

