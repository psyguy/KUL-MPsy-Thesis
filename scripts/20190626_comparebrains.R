rm(list = ls())

# source("./functions/functions_trial.R")
source("./functions/functions_netmeas.R")
source("./functions/functions_partition.R")

# for loop over all brain cases -------------------------------------------

sampled.path <- "./data/"
sampled.names <- list.files(path = sampled.path, pattern = "*.RData")
r.this <- sampled.names[grepl("_r-1_",sampled.names)]

# coefs.wb.b <- NULL
t <- Sys.time()
for(sampled in r.this){
  rm(brain_case)
  load(paste(sampled.path, sampled, sep = ""))
  brain_case@name %>% print()
  
  t1 <- Sys.time()
  for(i in 1:(brain_case@age$rewires-1)){
    m <- brain_case@history$mat.connectivity[[i]]
    if(is.null(m)) next
    paste("Adding coefs of rewiring", i) %>% print()
    coefs.wb.b <- coefs.wb.b %>% rbind(m %>% netmeas_wbcoefs(parameters = brain_case@parameters,
                                                             name = brain_case@name,
                                                             rewiring = i)
                                       )
  }
  Sys.time() - t1
}
Sys.time()-t

# 
# denom <- c(300, 50, 250)
# denom <- (denom*(denom-1)/2) %>% c(50*250)
# coefs.wb.b$Degree <- coefs.wb.b$Degree / rep(denom, nrow(coefs.wb.b)/4)
# 
# colnames(coefs.wb.b)[13] <- "Edge Density"

save_vars("coefs.wb.b", prefix = "coefs.wb.r1.5200")

# plotting coefficients over time -----------------------------------------
save_vars(list.of.vars = "coefs_all", prefix = "hpcJune26Harvest_tmp")
coefficient.name <- (coefs_all %>% colnames())[-1:-6]

# c <- coefs_all %>% filter(rewiring==1000) #filter(alphabeta.eps=="(0.5, 0.5)")

coefs_all_mutated <- coefs_all %>% 
  mutate(alphbet = paste0(alphabeta.eps,alphabeta.a))#,round)) %>% 
  filter(round==2)

for(i in coefficient.name){
  this.coefficient <- coefs_all_mutated %>%
    # filter(name == "Noa Adam") %>% 
    # filter(alphabeta.eps=="(5, 1)", alphabeta.a=="(5, 1)") %>% 
    select(i, rewiring, alphabeta.a, alphabeta.eps, name, alphbet)
  colnames(this.coefficient)[1] <- "value"
  # this.coefficient$value <- this.coefficient$value/co$value[1]
  coef.plot <- this.coefficient %>% ggplot(aes(x = rewiring,
                                               y = value,
                                               colour = alphbet#,
                                               # linetype = alphabeta.eps
                                               )) +
    geom_line(size = .75, alpha = 0.7) + 
    ggplot2::ylim(0,NA) +
    ggtitle(i) + theme(legend.position ="bottom")# c(0.7, 0.5))
    # guides(fill=guide_legend(nrow=2,byrow=TRUE))
  coef.plot %>% print()
  # paste0(i, "_uniques.png") %>% ggsave()
}



diamonds_reshaped <- data.frame(price = diamonds$price,
                                independent.variable = c(diamonds$carat,diamonds$cut,diamonds$color,diamonds$depth),
                                Clarity = rep(diamonds$clarity,times=4),
                                Variable.name = rep(c("Carat","Cut","Color","Depth"),each=nrow(diamonds)))

ggplot(diamonds_reshaped,aes(independent.variable,price,colour=Clarity)) + 
  geom_point(size=2) + facet_wrap(~Variable.name,scales="free_x") + 
  xlab("")


# Making coefs_all and plot pretty ----------------------------------------

coefs <- coefs_all_mutated #%>% colnames()
colnames(coefs) <- c("Owner", "Seed",
                     "Round", "Epsilon Proportion",
                     "a Proportion", "Rewiring",
                     "Clustering", "Efficiency",
                     "Small World", "Modularity",
                     "Avg Path Length", "Proportions(eps)(a)")

# coefs.long <- coefs %>% gather()

library(ggpubr)



for(i in 1:10){
  coefs.this.round <- coefs.wb.b %>%
    filter(Round==3) %>%
    filter(`.id` == "whole") %>% 
    filter(Rewiring>8500) %>% 
    filter(Rewiring<10000)
  # coefs.this.round <- coefs.this.round[5000:12500,]
  owner.name <- coefs.this.round$Owner[1] %>% as.character()
  title <- paste0("Coefficients for round ",
                  i,
                  " (",
                  owner.name,
                  ")")
  p1 <- ggplot(data = coefs.this.round,
               aes(x = Rewiring,
                   y = Clustering,
                   colour = Owner)) +
    geom_line(size = .75, alpha = 0.7) +
    ggplot2::ylim(0,NA)
  # ggplot2::ylim(0,NA)
  # p1
  
  p2 <- ggplot(data = coefs.this.round,
               aes(x = Rewiring,
                   y = `Small World`,
                   colour = Owner)) +
    geom_line(size = .75, alpha = 0.7) +
    ggplot2::ylim(0,NA)
  # p2
  
  p3 <- ggplot(data = coefs.this.round,
               aes(x = Rewiring,
                   y = Modularity,
                   colour = Owner)) +
    geom_line(size = .75, alpha = 0.7) + 
    ggplot2::ylim(0,NA)
  # p3
  
  p4 <- ggplot(data = coefs.this.round,
               aes(x = Rewiring,
                   y = `Avg Path Length`,
                   colour = Owner)) +
    geom_line(size = .75, alpha = 0.7) + 
    ggplot2::ylim(0,NA)
  
  figure <- ggarrange(p1, p2, p3, p4,
                      ncol=2, nrow=2,
                      common.legend = TRUE,
                      legend="bottom"
  )
  
  annotate_figure(figure, top = title)# %>% print()
  # paste0(title, ".png") %>% ggsave(width = 9,
  #                                  dpi = "retina")
}



# looking at connectivities -----------------------------------------------

library(seriation)

load("I:/Thesis_Codes/Thesis_R_Codes/data/eps-0v6v0_a-0v5v1_r-3_g-0.3k-5.2k_David De Ridder_20190626_1905.RData")
r.this.olive <- brain_case
load("I:/Thesis_Codes/Thesis_R_Codes/data/eps-1v5v0_a-0v6v0_r-3_g-0.3k-5.2k_Chiara Bosmans_20190626_1908.RData")
r.this.pink <- brain_case
load("I:/Thesis_Codes/Thesis_R_Codes/data/eps-0v5v1_a-0v6v0_r-3_g-0.3k-5.2k_Daan Pauwels_20190626_1904.RData")
r.this.red <- brain_case
rm(brain_case)

m <- r.this.pink@history$mat.connectivity[10400][[1]]
m[1:50,1:50] <- m[1:50,1:50]*3
m[1:50,51:300] <- m[1:50,51:300]*2
m[51:300,1:50] <- m[51:300,1:50]*2

# https://rstudio-pubs-static.s3.amazonaws.com/3486_79191ad32cf74955b4502b8530aad627.html
title <- paste0(b@name, " at ", b@age$rewires/1000, "k")
pimage(m,
       col = c("white", "brown3",
               "darkorchid1", "blue2"),
       key = FALSE,
       main = paste(title, "(unserialized)"))

pimage(m,
       seriate(m),
       col = c("white", "brown3",
               "darkorchid1", "blue2"),
       key = FALSE,
       main = paste(title, "(serialized)"))

size.minority <- (nrow(m)/6) %>% round(0)
size.majority <- (nrow(m)*5/6) %>% round(0)

g <- m %>% graph_from_adjacency_matrix(mode="undirected")
V(g)$partition <- c(rep("minority", size.minority),
                    rep("majority", size.majority))
# net <- g
# colrs <- c("gray50", "tomato", "gold")
# 
# V(net)$color <- colrs[V(net)$partition=="majority"]

# ceb <- cluster_edge_betweenness(net) 
# l <- g %>% layout_with_fr()
pg.9400 <- g %>% plot(#ceb,
           vertex.size=4,
           vertex.color=c("skyblue", "pink")[1+(V(g)$partition=="majority")],
           # layout = l,
           vertex.label = NA)


# coloring the edges ------------------------------------------------------
# c("white", "red", "black","blue") # 
# rain <- rainbow(4, alpha=.5)
# V(g)$color <- rain[V(g)$partition]
# 
# E(g)$color <- apply(as.data.frame(get.edgelist(g)), 1, 
#                     function(x){
#                       o <- rain[4]
#                       if(V(g)$partition[x[1]] == "minority" & V(g)$partition[x[2]] == "minority") o <-  rain[2]
#                       if(V(g)$partition[x[1]] == "majority" & V(g)$partition[x[2]] == "majority") o <- rain[3]
#                       o %>% return()
#                     })
# plot(g, vertex.size=4, vertex.label=NA, edge.color=E(g)$color)
