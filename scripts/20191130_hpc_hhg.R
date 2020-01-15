#!/usr/bin/env Rscript

Args <- commandArgs(TRUE)
counter <- as.numeric(as.character(Args[1]))


# rm(list=ls())
#source("./functions/functions_extract.R")
load("./data/StatusQuo_20191129_0425.RData")

#### NetSimile #################

library(tidyverse)
library(dplyr)
library(plyr)

library(HHG)
# writing nested for loop of HHG ------------------------------------------

d <- distances.connectivities

i.and.j <- list()
k <- 1
for(i in 1:length(d)){
  for(j in (i):length(d)){
    i.and.j[[k]] <- c(i,j)
    k <- k + 1
  }
}

i <- i.and.j[[counter]][1]
j <- i.and.j[[counter]][2]

hhg.results <- list(indices = c(i,j),
                    names = c(names(d[i]),names(d[j])),
                    hhg.values = hhg.test(d[[i]], d[[j]], nr.perm = 20000)
                    )

save_vars("hhg.results",
          prefix = paste("HHG", i, j, sep = "-")
          )

