## plotting adj+net with "PCA_angle" seriation method. Jesus.

source("functions/functions_my.R")
source("./functions/functions_extract.R")
load("./data-pc/snp_only-1e+6_20190723_1404.RData")

families <<- c("OC", "OT", "BL", "UT", "UC") %>% rep(each = 10)

make.df <- function(input.df, col){
  
  d <- input.df %>% 
    select("Partition",
           "Owner",
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
                           "Owner",
                           "Verbal.Description",
                           col,
                           "Club Size")
  output.df %>% return()
}

snp.new <- snp %>% 
  filter(Partition == "whole") %>% 
  mutate(fn = paste(Verbal.Description, Owner, sep="_")) %>%
  arrange(fn) %>%
  mutate(fn = paste0(rep(c(3,2,1,4,5),each=10), "_",fn)) %>% 
  arrange(fn) %>% 
  cbind(family.codes = families, num =1:50) %>% 
  mutate(newcodes = paste0(num,"_",family.codes,rep(c(1:10), 5), "_", Owner)) %>% 
  mutate(famcode.num = paste0(family.codes, rep(c(1:10), 5)))



all.owners <- snp.new$famcode.num %>% 
  unique() %>% 
  sort() %>% 
  as.character() %>%
  as.list()



for(i in 1:nrow(snp.new)){
  m <- snp.new$adj.mat.vect[[i]] %>% vec2mat()
  title <- snp.new$famcode.num[[i]]
  print(paste (i, title))
  
  extract_plotcon(m, title)
  extract_plotnet(m, title)
  
}

