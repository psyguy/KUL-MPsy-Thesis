# What has been so far used (dis)similarity with sloppiness. Here I am correcting the
# definitions and cleaning the code. Moreover, I will plot resemblences/differentiations
# as networks
rm(list=ls())
load("./data/signature-and-HHG_20191202_1704.RData")
load("./data/HHG-unscaled_20191216_1510.RData")
source("./scripts/20191202_family-comparison-heatmaps.R")

library(plyr)
library(viridis)
library(igraph)
col.pallett <- inferno(10,1,0,1,1)


# functions used for plotting heatmaps ------------------------------------

sim.hm <- function(m, title = "", col = col.pallett){
  f <- families
  m <- m %>% as.matrix() %>% unname()
  Heatmap(m,
          col = col,
          column_title = title,
          name = "Similarities",
          na_col = "black",
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          row_split = factor(f, levels = unique(f)),
          column_split = factor(f, levels = unique(f)),
          row_gap = unit(0.5, "mm"), column_gap = unit(0.5, "mm"), border = TRUE)
  }

f.corrplot <- function(m,
                       type = "lower",
                       col = col.pallett) corrplot(m,
                                                   method = "shade",
                                                   type = type,
                                                   is.corr = T,
                                                   col = rep(col,2),
                                                   addCoef.col = col.pallett[1],
                                                   diag = TRUE,
                                                   cl.pos = "n",
                                                   tl.pos = "n",
                                                   number.cex = .6)

# making Sim/Diff matrices ------------------------------------------------

hhg.act <- df.hhg.act %>% make.hhg.mat()
hhg.con <- df.hhg.con %>% make.hhg.mat()
l.Similarity <- list(`Anatomical HHG similarity` = 1 - hhg.con,
                     `Anatomical NetSimile similarity` = 1 - make.halved.mats(s.con, s.con, TRUE),
                     `Functional HHG similarity` = 1 - hhg.act,
                     `Functional NetSimile similarity` = 1 - make.halved.mats(s.act, s.act, TRUE)
                     )

l.Similarities.halved <- list(`Anatomical` = make.hhg.s(l.Similarity$`Anatomical HHG similarity`,
                                                          l.Similarity$`Anatomical NetSimile similarity`),
                              `Functional` = make.hhg.s(l.Similarity$`Functional HHG similarity`,
                                                      l.Similarity$`Functional NetSimile similarity`)
                              )
l.Resemblance <- l.Similarity %>%
  map(~make.block.mean(make.halved.mats(.,.,TRUE), TRUE))

for(l in l.Resemblance) rownames(l) <- colnames(l) <- unique(families)

Differentiation <- l.Resemblance %>% 
  ldply(function(l){
    sm <- NULL
    for(i in 1:5){
      # sm <- 
      x <- l[i,]
      sm[i] <- 1 - x[i]/mean(x[-i])
    }
    names(sm) <- families %>% unique()
    sm %>% return()
  }
  )

g.d <- Differentiation %>% gather("Family", "Differentiation", -.id)

g.d$Family <- factor(g.d$Family, levels = unique(g.d$Family))

colnames(g.d)[1] <- "Resemblance measure"

l.Similarities.halved$Anatomical %>% sim.hm()
l.Similarities.halved$Functional %>% sim.hm()

l.Resemblance$`Anatomical HHG similarity` %>% f.corrplot()
l.Resemblance$`Anatomical NetSimile similarity` %>% f.corrplot(type = "upper")
l.Resemblance$`Functional HHG similarity` %>% f.corrplot()
l.Resemblance$`Functional NetSimile similarity` %>% f.corrplot(type = "upper")

# Plotting Similarity/Resemblance/Differentiation -------------------------

g.d %>% ggplot(aes(x = `Resemblance measure`,
                   y = Differentiation,
                   fill = Family),
               title = "avg") +
  ylim(-.5,.5) +
  geom_bar(position = "dodge", size = 1, width = 0.5, stat="identity")


# Plotting graphs ---------------------------------------------------------

## Change the index in this line
plot.name <- pull(Differentiation, .id)[4]

differ <- Differentiation %>%
  filter(.id == plot.name) %>% 
  select(-.id) %>% 
  as.numeric()

positive.differ <- rep("*", 5)
positive.differ[differ<0] <- ""

resemb <- l.Resemblance[[1]]
colnames(resemb) <- families %>%
  unique() %>% 
  paste0(positive.differ,.)

g <- resemb %>% graph_from_adjacency_matrix(mode = "undirected",
                                            weighted = TRUE,
                                            diag = FALSE)

#Color scaling function
c_scale <- colorRamp(col.pallett)

#Applying the color scale to edge weights.
#rgb method is to convert colors to a character vector.
E(g)$width <- E(g)$weight
E(g)$color <- apply(c_scale(E(g)$weight), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
V(g)$color <- apply(c_scale(diag(resemb)), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
# V(g)$size <- 100*abs(differ)

v.frame.color <- rep("green", 5)
v.frame.color[differ<0] <- "white"


radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}

la <- layout.circle(g)

par(mar=c(0,0,0,0)+1.7)
g %>% plot(edge.width = 30*E(g)$weight,
           main = gsub(" similarity", "", plot.name),
           vertex.size = 300*abs(differ),
           vertex.label.dist = 3.0,
           # vertex.frame.color = v.frame.color,
           vertex.frame.width = 100,
           vertex.label.degree = radian.rescale(x=1:5, direction=-1, start=6),
           layout=la
           )
