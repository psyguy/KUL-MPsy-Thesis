# Notes -------------------------------------------------------------------
##########################################################
######## plotting adj+net with PCA_angle seriation #######
######## 2021-06-14 ######################################
##########################################################
## 
## From commit 2aa903a on, the following snp data files should be used:
##    - snp_all_5k_20210628_1238.RData (zipped into three parts, for GitHub)
##    - snp.lean_all_5k_20210628_1241.RData
##    - snp_only-1e+6_20210628_1241.RData
## 
##########################################################


# Loading functions and snp -----------------------------------------------

source("functions/functions_my.R")
source("./functions/functions_extract.R")
load("./data/snp_only-1e+6_20210628_1241.RData")


make.df <- function(input.df, col){
  
  d <- input.df %>% 
    select("Partition",
           "ID",
           "Verbal.Description",
           col)
  
  output.df <- NULL
  for(p in 1:4){
    this.df <- select(d, -col)[p,]
    vec <- pull(d, col)[[p]]
    len <- length(vec)
    x <- this.df[rep(seq_len(nrow(this.df)),
                     each = len),]
    
    output.df <- output.df %>%
      rbind(
        cbind(x, vec, c(1:len))
      )
  }
  
  colnames(output.df) <- c("Partition",
                           "ID",
                           "Verbal.Description",
                           col,
                           "Club Size")
  output.df %>% return()
}

snp.new <- snp %>% 
  filter(Partition == "whole")
  # The rest is redundant, as they have been corrected in 2aa903a
  # %>% 
  # mutate(fn = paste(Verbal.Description, Owner, sep="_")) %>%
  # arrange(fn) %>%
  # mutate(fn = paste0(rep(c(3,2,1,4,5),each=10), "_",fn)) %>% 
  # arrange(fn) %>% 
  # cbind(family.codes = families, num =1:50) %>% 
  # mutate(newcodes = paste0(num,"_",family.codes,rep(c(1:10), 5), "_", Owner)) %>% 
  # mutate(famcode.num = paste0(family.codes, rep(c(1:10), 5)))

all.owners <- snp.new$fn %>% 
  unique() %>% 
  sort() %>% 
  as.character() %>%
  as.list()

i <- 2

for(i in 1:nrow(snp.new)){
  m <- snp.new$adj.mat.vect[[i]] %>% vec2mat()
  title <- snp.new$ID[[i]]
  print(paste (i, title))
  
  plotcon <- extract_plotcon(m, title)
  plotnet <- extract_plotnet(m, title)
  
}

librarian::shelf(kwb-r/kwb.plot)
