#!/usr/bin/env Rscript

# Args <- commandArgs(TRUE)
# rounds <- as.numeric(Args[1])
rm(list = ls())

source("./functions/functions_extract.R")

pattern <- "snapshots_25e3_"

paste("Now reading round", r_, "of", pat.tmp, "index", rounds) %>% print()

path.to.snapshots <- "./data/snapshots/"

sampled.names <- list.files(path = path.to.snapshots, pattern = "*.RData")
r.this <- sampled.names[grepl(pattern, sampled.names)]
  # r.this <- r.this[grepl("_g-", r.this)]
  # r.this <- r.this[!grepl("coefs.wb.", r.this)]
  
snapshots_all_25e3 <- NULL
t <- Sys.time()
for(sampled in r.this){
  load(paste(path.to.snapshots, sampled, sep = ""))
  snapshots_all_25e3 <- snapshots_all_25e3 %>%
    rbind(snapshots_25e3)
}
(Sys.time()-t) %>% paste("for", pattern) %>% print()

snp <- snapshots_all_25e3

coefs.names <- colnames(snp)[9:15]

color.bg <- "white"
color.mino <- "deepskyblue3" #"blue2"
color.majo <- "orangered" #"brown3"
color.inter <- "olivedrab2" #"green4"
color.whole <- "azure4"

max.cl <- snp$Clustering %>% max()
max.sw <- snp$Small.World %>% max()
max.mo <- snp$Modularity %>% max()
max.as <- snp$Assortativity %>% max()
min.as <- snp$Assortativity %>% min()
max.rc <- snp$Rich.Club %>% max()
# max.de <- snp$Degree %>% max()
max.pl <- snp$Avg.Path.Length %>% max()
max.ef <- snp$Efficiency %>% max()

name <- "Logan Wauters"#title# <- snp$Owner %>% as.character() %>% unique()
t3 <- Sys.time()
for(name in names){
  print(name)
  coefs.this.round <- snp %>% filter(Owner == name)
  owner.name <- coefs.this.round$Owner[1] %>% as.character()
  title <- paste0("Evolution of coefficients for ",
                  name,
                  "eps ",
                  coefs.this.round$`Epsilon Proportion`,
                  ", a ",
                  coefs.this.round$`a Proportion`)
  
  p1 <- ggplot(data = coefs.this.round,
               aes(x = Rewiring,
                   y = Clustering,
                   colour = Partition)) +
    geom_line(size = .75, alpha = 0.7) +
    scale_colour_manual(values = c(color.inter, color.majo,
                                   color.mino, color.whole)) +
    ggplot2::ylim(0,max.cl)
  
  p2 <- ggplot(data = coefs.this.round,
               aes(x = Rewiring,
                   y = Small.World,
                   colour = Partition)) +
    geom_line(size = .75, alpha = 0.7) +
    scale_colour_manual(values = c(color.inter, color.majo,
                                   color.mino, color.whole)) +
    ggplot2::ylim(0,max.sw)
  
  p3 <- ggplot(data = coefs.this.round,
               aes(x = Rewiring,
                   y = Modularity,
                   colour = Partition)) +
    geom_line(size = .75, alpha = 0.7) + 
    scale_colour_manual(values = c(color.inter, color.majo,
                                   color.mino, color.whole)) +
    ggplot2::ylim(0,max.mo)
  
  p4 <- ggplot(data = coefs.this.round,
               aes(x = Rewiring,
                   y = Assortativity,
                   colour = Partition)) +
    geom_line(size = .75, alpha = 0.7) + 
    scale_colour_manual(values = c(color.inter, color.majo,
                                   color.mino, color.whole)) +
    ggplot2::ylim(min.as,max.as)
  
  p5 <- ggplot(data = coefs.this.round,
               aes(x = Rewiring,
                   y = Avg.Path.Length,
                   colour = Partition)) +
    geom_line(size = .75, alpha = 0.7) + 
    scale_colour_manual(values = c(color.inter, color.majo,
                                   color.mino, color.whole)) +
    ggplot2::ylim(0,max.pl)
  
  p6 <- ggplot(data = coefs.this.round,
               aes(x = Rewiring,
                   y = Rich.Club,
                   colour = Partition)) +
    geom_line(size = .75, alpha = 0.7) + 
    scale_colour_manual(values = c(color.inter, color.majo,
                                   color.mino, color.whole)) +
    ggplot2::ylim(0,max.ef)
  
  
  figure <- ggarrange(p1, p2, p3, p4, p5, p6,
                      ncol=1, #nrow=2,
                      common.legend = TRUE,
                      legend="bottom"
  )
  
  annotate_figure(figure, top = title)# %>% print()
  paste0("wb_", name, ".png") %>% ggsave(width = 7,
                                         height = 42,
                                         dpi = "retina")
}
Sys.time() - t3


# plots graphs with colored edges -----------------------------------------


sampled.path <- "./data/"
sampled.names <- list.files(path = sampled.path, pattern = "*.RData")
r.this <- sampled.names[grepl("Logan",sampled.names)]
load(paste(sampled.path, r.this, sep = ""))

rew <- 100000

t4 <- Sys.time()
m <- brain_case@history$mat.connectivity[[rew]]
g <- m %>% graph_from_adjacency_matrix(mode="undirected")

title <- brain_case@name %>% paste("after", rew, "rewirings")
# num.mino <- (brain_case@parameters$n_edges/6) %>% round(digits=0)

V(g)$partition <- c(rep("minority", 50),
                    rep("majority", 250))

# select edges and set color 
E(g)[V(g)[partition == "minority"] %--% V(g)[partition == "majority"]]$color <- color.inter
E(g)[V(g)[partition == "majority"] %--% V(g)[partition == "majority"]]$color <- color.majo
E(g)[V(g)[partition == "minority"] %--% V(g)[partition == "minority"]]$color <- color.mino
# plot
set.seed(1)
paste0(title, "_network", ".png") %>%
  png(width = width.column.report*2, height = width.column.report*2, res = 200)
#Set the margin size (huge margins)
par(mar=c(.1,.1,.1,.1))
g %>% plot(vertex.size = 0,
           # vertex.color = c("skyblue", "red")[1 + (V(g)$partition == "majority")],
           add=TRUE,
           vertex.label = NA,
           vertex.color = "black",
           edge.width = 1,
           # rescale = F,
           # xlim = c(-5,5),
           edge.curved= 0.5#,
           # main = paste(brain_case@name,
           #              "at",
           #              rew,
           #              "rewirings")
)
dev.off()


s <- m %>% seriate()
m[1:50,1:50] <- m[1:50,1:50]*3
m[1:50,51:300] <- m[1:50,51:300]*2
m[51:300,1:50] <- m[51:300,1:50]*2

width.column.report <- 2240
col.pimage <- c(color.bg, color.majo,
                color.inter, color.mino)


paste0(title, "_unserialized", ".png") %>%
  pdf(width = width.column.report,
      height = width.column.report)#,
      res = 400)
pimage(m,
       col = col.pimage,
       key = FALSE#,
       # main = paste(brain_case@name,
       #              "at",
       #              rew,
       #              "rewirings")# (unserialized)")
)
dev.off()


paste0(title, "_serialized", ".png") %>%
  png(width = width.column.report, height = width.column.report, res = 400)
pimage(m,
       s,
       col = col.pimage,

       key = FALSE#,
       # main = paste(brain_case@name,
       #              "at",
       #              rew,
       #              "rewirings")# (serialized)")
)
dev.off()

# putting adj mats next to each other
gl <- rep(title,2) %>%
  paste0(c("_unserialized.png", "_serialized.png")) %>% 
  as.list() %>%
  lapply(png::readPNG) %>% 
  lapply(grid::rasterGrob)

paste0(title, "_connectivities", ".png") %>%
  png(width = width.column.report*2, height = width.column.report, res = 400)
gr.con <- gridExtra::grid.arrange(grobs=gl, ncol=2,
                                  # right = paste(brain_case@name,
                                  #               "at",
                                  #               rew,
                                  #               "rewirings"),
                                  padding = unit(1, "point")
)
dev.off()


conns_and_net <- rep(title,2) %>%
  paste0(c("_connectivities.png", "_network.png")) %>% 
  as.list() %>%
  lapply(png::readPNG) %>% 
  lapply(grid::rasterGrob)
require(grid)
paste0(title, "", ".png") %>%
  png(width = width.column.report*2, height = width.column.report*3, res = 400)
gr.net <- gridExtra::grid.arrange(grobs=conns_and_net, nrow=2,
                                  top = textGrob(paste("\n",
                                                       brain_case@name,
                                                       "after",
                                                       rew,
                                                       "rewirings"),
                                                 gp=gpar(fontsize=40,font=8)
                                  )
)
dev.off()
Sys.time()-t4




