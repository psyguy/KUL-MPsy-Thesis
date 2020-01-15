# rm(list=ls())

library(tidyverse)

hhg.files <- list.files(path = "./data/HHG/")
cn <- c(
  "index.1",
  "index.2",
  "name.1",
  "name.2",
  "sum.chisq",
  "sum.lr",
  "max.chisq",
  "max.lr",
  "perm.pval.hhg.sc",
  "perm.pval.hhg.sl",
  "perm.pval.hhg.mc",
  "perm.pval.hhg.ml",
  "stat.type",
  "n")
df.hhg.con <- df.hhg.act <- data.frame(matrix(ncol = length(cn), nrow = 0))
for(i in hhg.files){
  load(paste0("./data/HHG/",i))
  df.hhg.con <- hhg.con %>% unlist() %>% as.character() %>% rbind(df.hhg.con,.) %>% mutate_all(as.character)
  df.hhg.act <- hhg.act %>% unlist() %>% as.character() %>% rbind(df.hhg.act,.) %>% mutate_all(as.character)
}

colnames(df.hhg.con) <- colnames(df.hhg.act) <- cn

df.h <- df.hhg.con %>% select(index.1:index.2,sum.chisq:perm.pval.hhg.ml) %>% mutate_all(as.numeric)
h.mat <- matrix(nrow = 50, ncol = 50)
for(i in 1:nrow(df.h)){
  h <- df.h[i,]
  h.mat[h$index.1, h$index.2] <- df.hhg.act[i,]$perm.pval.hhg.ml %>% as.numeric() #1 - log(h$sum.chisq)/max(log(df.h$sum.chisq))
  h.mat[h$index.2, h$index.1] <- h$perm.pval.hhg.ml#%>% log()
  }
# diag(h.mat) <- NaN
(h <- h.mat) %>% pimage


# colors = list(
#   bg = "white",
#   mino = "deepskyblue3", #"blue2"
#   majo = "orangered", #"brown3"
#   inter = "olivedrab2", #"green4"
#   whole = "azure4"
# )

# plotting with heatmap ---------------------------------------------------

# library(devtools)
# install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
diag(h) <- 1

families <- df.hhg.act %>%
  mutate(index = as.numeric(index.1)) %>% 
  arrange(index) %>%
  pull(name.1) %>%
  unique() %>% 
  gsub(" .*", "", .)

library(circlize)
Heatmap(s,
  # col = colorRamp2(c(0, 0.05, 1),c("orangered", "azure4", "white")),
  col = c("white", "deepskyblue3"),
  na_col = "black",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = families,
  column_split = families,
  row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border = TRUE)



# doing similar things with NetSimile signature distances -----------------

# Euclidean distances
s.euc.con <- signatures.connectivities[-1] %>%
  dist() %>% 
  as.matrix()
s.euc.act <- signatures.activities[-1] %>%
  dist() %>% 
  as.matrix()
# Canberra distances
s.con <- signatures.connectivities[-1] %>%
  dist(method = "canberra") %>% 
  as.matrix()
s.act <- signatures.activities[-1] %>%
  dist(method = "canberra") %>% 
  as.matrix()
s.con[upper.tri(s.con)] <- 0
s.act[lower.tri(s.act)] <- 0
s.act[is.na(s.act)] <- 0


s <- s.con/max(s.con) + s.act/max(s.act)
s[is.na(s)] <- 0
s[s>1] <- 1



h.c <- h
h.c[upper.tri(h.c)] <- 0
h.c <- h.c + t(h.c)
diag(h.c) <- NA
mean.h.c <- matrix(nrow = 50, ncol = 50) 

for(i in 1:5){
  for(j in 1:5){
    mean.h.c[(10*i-9):(i*10),(10*j-9):(j*10)] <- mean(h.c[(10*i-9):(i*10),(10*j-9):(j*10)], na.rm = TRUE)
  }
}

s.con <- signatures.connectivities
s.act <- signatures.connectivities

save_vars(c("s.con", "s.euc.con", "s.act", "s.euc.act", "df.hhg.act", "df.hhg.con"), prefix = "signature-and-HHG")
# save_vars(c("df.hhg.act", "df.hhg.con"), prefix = "HHG-unscaled")
