rm(list = ls())

# source("./functions/functions_trial.R")
source("./functions/functions_netmeas.R")
source("./functions/functions_partition.R")

load("I:/Thesis_Codes/Thesis_R_Codes/data/eps-0v6v0_a-0v5v1_r-3_g-0.3k-5.2k_David De Ridder_20190626_1905.RData")
r3.olive <- brain_case
load("I:/Thesis_Codes/Thesis_R_Codes/data/eps-1v5v0_a-0v6v0_r-3_g-0.3k-5.2k_Chiara Bosmans_20190626_1908.RData")
r3.pink <- brain_case
load("I:/Thesis_Codes/Thesis_R_Codes/data/eps-0v5v1_a-0v6v0_r-3_g-0.3k-5.2k_Daan Pauwels_20190626_1904.RData")
r3.red <- brain_case
rm(brain_case)

# trying within/between partition measures --------------------------------

b <- r3.olive
# b <- r3.pink
# b <- r3.red
p <- b@parameters
name <- b@name

# coefs.wb.b <- NULL
t1 <- Sys.time()
for(i in 1:b@age$rewires){
  m <- b@history$mat.connectivity[[i]]
  if(is.null(m)) next
  paste("Adding coefs of rewiring", i) %>% print()
  coefs.wb.b <- coefs.wb.b %>% rbind(m %>% netmeas_wbcoefs(parameters = b@parameters,
                                                           name = b@name,
                                                           rewiring = i)
                                     )
}
Sys.time() - t1

save_vars("coefs.wb.b", prefix = "round3 partitioned measures")

partition.list <- c("minority", "majority",
                    "interpartition", "whole")

max.cl <- coefs.wb.b$Clustering %>% max()
max.sw <- coefs.wb.b$`Small World` %>% max()
max.mo <- coefs.wb.b$Modularity %>% max()
max.pl <- coefs.wb.b$`Avg Path Length` %>% max()

for(partition in partition.list){
  print(partition)
  coefs.this.round <- coefs.wb.b %>% filter(`.id` == partition)
  owner.name <- coefs.this.round$Owner[1] %>% as.character()
  title <- paste0("Coefficients for ",
                  partition)#,
                  # " (",
                  # owner.name,
                  # ")")
  p1 <- ggplot(data = coefs.this.round,
               aes(x = Rewiring,
                   y = Clustering,
                   colour = Owner)) +
    geom_line(size = .75, alpha = 0.7) +
    ggplot2::ylim(0,max.cl)
    # ggplot2::ylim(0,NA)
  # p1
  
  p2 <- ggplot(data = coefs.this.round,
               aes(x = Rewiring,
                   y = `Small World`,
                   colour = Owner)) +
    geom_line(size = .75, alpha = 0.7) +
    ggplot2::ylim(0,max.sw)
  # p2
  
  p3 <- ggplot(data = coefs.this.round,
               aes(x = Rewiring,
                   y = Modularity,
                   colour = Owner)) +
    geom_line(size = .75, alpha = 0.7) + 
    ggplot2::ylim(0,max.mo)
  # p3
  
  p4 <- ggplot(data = coefs.this.round,
               aes(x = Rewiring,
                   y = `Avg Path Length`,
                   colour = Owner)) +
    geom_line(size = .75, alpha = 0.7) + 
    ggplot2::ylim(0,max.pl)
  
  figure <- ggarrange(p1, p2, p3, p4,
                      ncol=2, nrow=2,
                      common.legend = TRUE,
                      legend="bottom"
  )
  
  annotate_figure(figure, top = title)# %>% print()
  paste0(title, ".png") %>% ggsave(width = 9,
                                   dpi = "retina")
}


m <- r3.pink@history$mat.connectivity[rew][[1]]
m[1:50,1:50] <- m[1:50,1:50]*3
m[1:50,51:300] <- m[1:50,51:300]*2
m[51:300,1:50] <- m[51:300,1:50]*2

# https://rstudio-pubs-static.s3.amazonaws.com/3486_79191ad32cf74955b4502b8530aad627.html
title <- paste0(b@name, " at ", b@age$rewires)

paste0(title, "_1", ".png") %>% 
  png(width = 1440, height = 1440, res = 200)
pimage(m,
       col = c("white", "brown3",
               "darkorchid1", "blue2"),
       key = FALSE)#,
       # main = paste(title, "(unserialized)"))
dev.off()


paste0(title, "_2", ".png") %>% 
  png(width = 1440, height = 1440, res = 200)
pimage(m,
       seriate(m),
       col = c("white", "brown3",
               "darkorchid1", "blue2"),
       key = FALSE)#,
       # main = paste(title, "(serialized)"))
dev.off()


rl = lapply(sprintf(paste0(title,"_%i.png"),
                    1:2), png::readPNG)
gl = lapply(rl, grid::rasterGrob)
gridExtra::grid.arrange(grobs=gl, ncol=2)

