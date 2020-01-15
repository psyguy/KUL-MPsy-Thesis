# rm(list=ls())

load("data-pc/snp_only-1e+6_20190723_1404.RData")
source("functions/functions_my.R")
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
  filter(!(num %in%c(12,13,41,43))) # removing the dead brains


tmp <- snp.new %>% filter(num == 1)
extract_RCnorm_sig <- function(tmp, num.rand = 2, sig.level = 0.05){
  
  max.club.size <- 200
  rc_norm <- tmp$adj.mat.vect %>% 
    vec2mat() %>% 
    graph_from_adjacency_matrix() %>% 
    rich_club_norm(num.rand)
  
  rc_norm <- rc_norm[1:max.club.size,]
  
  rc.res <- matrix(nrow = max.club.size, ncol = 5) %>% as.data.frame()
  colnames(rc.res) <- c("num", "family.codes", "Club Size", "Rich Club", "p.value")
  
  rc.res$num <- tmp$num %>% rep(max.club.size)
  rc.res$family.codes <- tmp$family.codes %>% rep(max.club.size)
  rc.res$`Club Size` <- 1:max.club.size
  rc.res$`Rich Club`[1:nrow(rc_norm)] <- rc_norm$norm
  rc.res$p.value[1:nrow(rc_norm)] <- rc_norm$p
  # rc.res$p.value[rc.res$p.value>sig.level] <- NA
  
  rc.res %>% return()
}

rc.res <- matrix(nrow = 0, ncol = 5) %>% as.data.frame()
colnames(rc.res) <- c("num", "family.codes", "Club Size", "Rich Club", "p.value")

for(i in 1:46){
  t1 <- Sys.time()
  xt <- snp.new[i,]
  rc.res <- rc.res %>% rbind(extract_RCnorm_sig(xt, 200))
  paste(i, "RC of", xt$newcodes, "took", Sys.time()-t1) %>% print()
}

save_vars("rc.res", prefix = "Normalized rich club of whole")



# plotting RC values and significance levels ------------------------------


rc.res$num <- rc.res$num %>% as.character() %>% as.factor()

fam.code <- "BL"

rc.this <- rc.res %>%
  mutate(rc.p = `Rich Club`) %>% 
  filter(family.codes== fam.code)

rc.this$rc.p[rc.this$p.value>0.05] <- NA
rc.this$`Rich Club`[rc.this$`Rich Club`==0] <- NA

rc.this %>%  ggplot(
            aes(x = `Club Size`,
                y = `Rich Club`,
                colour = num)) +
  geom_line(size = 1.5, alpha = 0.25) +
  ggplot2::xlim(0,150) +
  ggplot2::ylim(0.5,2.20) +
  ggtitle(paste("Normalized Rich Club for", fam.code)) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  geom_point(aes(x = `Club Size`,
                 y = rc.p,
                 colour = num),
             data = rc.this,
             size = .95,
             shape = 16) +
  geom_hline(yintercept = 1,linetype="dashed")



# Redoing family evolution plots, this time without RC/efficiency ---------

all.partitions <- snp$Partition %>% 
  unique() %>% 
  sort() %>% 
  as.character() %>%
  as.list()

all.vd <- snp$Verbal.Description %>% 
  unique() %>% 
  sort() %>% 
  as.character() %>%
  as.list()


system.time(
  for(vd in all.vd){
    for(part in all.partitions){
      extract_plotcoefs.glued(this.Verbal.Description = vd,
                              this.Partition = part, 
                              snp=snp,
                              coef.range = coefindex)
    }
  }
)

