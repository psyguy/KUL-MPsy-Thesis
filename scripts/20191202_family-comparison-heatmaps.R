rm(list=ls())

library(tidyverse)
library(ComplexHeatmap)
library(seriation)
library("corrplot")

load("./data/signature-and-HHG_20191202_1704.RData")
load("./data/HHG-unscaled_20191216_1510.RData")

families <<- df.hhg.act %>%
  mutate(index = as.numeric(index.1)) %>% 
  arrange(index) %>%
  pull(name.1) %>%
  unique() %>% 
  gsub(" .*", "", .)


# functions ---------------------------------------------------------------

make.hhg.mat <- function(df, hhg.value = "perm.pval.hhg.ml"){
  df.h <- df %>% mutate_at(c("index.1","index.2"),as.numeric)
  h.mat <- matrix(nrow = 50, ncol = 50)
  for(i in 1:nrow(df.h)){
    h <- df.h[i,]
    h.mat[h$index.1, h$index.2] <- h.mat[h$index.2, h$index.1] <- h %>%
      select(hhg.value) %>% as.numeric()
  }
  h.mat %>% return()
}

make.halved.mats <- function(m.lower, m.upper, standardize = FALSE){
  m.lower[upper.tri(m.lower)] <- 0
  m.upper[lower.tri(m.upper)] <- 0
  o <- m.lower + m.upper
  if(standardize) o <- m.lower/max(m.lower, na.rm = TRUE) + m.upper/max(m.upper, na.rm = TRUE)
  diag(o) <- NA
  o %>% return()
}

make.block.mean <- function(m, single.cell = FALSE){
  o <- matrix(nrow = 50, ncol = 50) 
  for(i in 1:5){
    for(j in 1:5){
      o[(10*i-9):(i*10),(10*j-9):(j*10)] <- mean(m[(10*i-9):(i*10),(10*j-9):(j*10)], na.rm = TRUE)
    }
  }
  if(single.cell){
    o <- matrix(nrow = 5, ncol = 5)
    for(i in 1:5){
      for(j in 1:5){
        o[i,j] <- mean(m[(10*i-9):(i*10),(10*j-9):(j*10)], na.rm = TRUE)
      }
    }
  }
  o %>% return()
}

make.together <- function(m){
  m <- make.halved.mats(m, m, TRUE)
  m %>% make.block.mean() %>% make.halved.mats(m,.) %>% return()
}


library(circlize)
hm <- function(m, title){
  f <- families
  if(nrow(m)==5) f <- families %>% unique()
  Heatmap(m,
          # col = colorRamp2(c(0, 0.05, 1),c("orangered", "azure4", "white")),
          col = c("white", "deepskyblue3"),
          row_title = title,
          na_col = "black",
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          row_split = factor(f, levels = unique(f)),
          column_split = factor(f, levels = unique(f)),
          row_gap = unit(0.5, "mm"), column_gap = unit(0.5, "mm"), border = TRUE)
}

make.hhg.s <- function(m.hhg, m.s) make.halved.mats(m.hhg, make.halved.mats(m.s, m.s, TRUE))

# making matrices ---------------------------------------------------------

# hhg.act <- df.hhg.act %>% make.hhg.mat()
# hhg.con <- df.hhg.con %>% make.hhg.mat()
# 
# l.family.sim <- list(`Connectivity HHG test` = hhg.con,
#                      `Connectivity NetSimile distances` = s.con,
#                      `Activity HHG test` = hhg.act,
#                      `Activity NetSimile distances` = s.act
#                      )
# 
# 
# (l.family.halved <- list(`Connectivity` = make.hhg.s(hhg.con,s.con),
#                         `Activity` = make.hhg.s(hhg.act,s.act)
#                         )
#   )
# 
# l.avg <- l.family.sim %>% map(~make.block.mean(make.halved.mats(.,.,TRUE), TRUE))
# 
# 
# # l.avg <- list(`Connectivity HHG test` = make.halved.mats(hhg.con, hhg.con, TRUE) %>% make.block.mean(TRUE),
# #               `Connectivity NetSimile distances` = make.halved.mats(s.con, s.con, TRUE) %>% make.block.mean(TRUE))
# 
# 
# # l.family.sim %>% 
# 
# 
# l <- l.avg[[3]]
# rownames(l) <- colnames(l) <- unique(families)
# 
# l %>% corrplot(method = "square",
#                type = "lower",
#                is.corr = T,
#                # addgrid.col = NA,
#                # addCoefasPercent = T,
#                addCoef.col = "black",
#                tl.pos = "r",
#                tl.col = "black",
#                number.cex = .6
#                )
# 
# # l %>% corrplot(method = "number",
# #                type = "upper",
# #                is.corr = T,
# #                # addgrid.col = NA,
# #                # addCoefasPercent = T,
# #                # addCoef.col = "black",
# #                tl.pos = "r",
# #                tl.col = "black",
# #                number.cex = .9,
# #                add = T
# #                )
# 
# 
# 
# 
# 
# 
# l[lower.tri(l, diag = T)]
# 
# corrplot()
# ## trial here, do not run
# ## .euc. distances yield less significant distinctions, hence use canberra
# # 
# # hhg.act %>% make.together()
# # s.euc.act %>% make.together()
# # s.act %>% make.together()
# # 
# # hhg.con %>% make.together() 
# # s.euc.con %>% make.together()
# # s.con %>% make.together()
# # 
# # 
# # s.euc.act %>% make.halved.mats(s.act) %>%  pimage
# # s.euc.con %>% make.halved.mats(s.con) %>%  pimage
# # s.euc.act %>% make.block.mean() %>% make.halved.mats(make.block.mean(s.act)) %>% pimage
# # s.euc.con %>% make.block.mean() %>% make.halved.mats(make.block.mean(s.con)) %>% pimage
# # hhg.con %>% make.block.mean() %>% make.halved.mats((hhg.act)) %>% pimage
# 




