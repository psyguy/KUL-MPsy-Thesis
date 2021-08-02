# What has been so far used (dis)similarity with sloppiness. Here I am correcting the
# definitions and cleaning the code. Moreover, I will plot resemblences/differentiations
# as networks

# *** NOTE *** th edissimilarities are being used in this file

# rm(list=ls())
load("./data/signature-and-HHG_20191202_1704.RData")
load("./data/HHG-unscaled_20191216_1510.RData")
source("./scripts/20191202_family-comparison-heatmaps.R")
load("./data/NetSimile-distances_20210802_1319.RData")

load("./data/snp_only-1e+6_20210628_1241.RData")
snp.micro <- snp %>% 
  filter(Partition == "whole") %>% 
  select(ID,
         Owner,
         Verbal.Description,
         fn) %>% 
  mutate(index = substring(fn, 1,2))

df.hhg.act <- df.hhg.act %>% 
  mutate(name.1 = gsub(".*_","", name.1),
         name.2 = gsub(".*_","", name.2),
         index.1 = gsub(".*_","", name.1),
         index.2 = gsub(".*_","", name.2))

df.hhg.con <- df.hhg.con %>% 
  mutate(name.1 = gsub(".*_","", name.1),
         name.2 = gsub(".*_","", name.2),
         index.1 = gsub(".*_","", name.1),
         index.2 = gsub(".*_","", name.2))


library(plyr)
library(viridis)
library(igraph)
col.pallett <- inferno(10,1,0,1,-1)

families <<- c("HC", "HD", "BL", "SD", "SC") %>% rep(each = 10)
# functions used for plotting heatmaps ------------------------------------

sim.hm <- function(m, title = "", col = col.pallett){
  f <- families
  m <- m %>% as.matrix()# %>% unname()
  colnames(m) <- rownames(m) <- rep(1:10,5) %>% paste0(families,.)
  Heatmap(m,
          col = col,
          column_title = title,
          name = "Dissimilarities",
          heatmap_legend_param = list(
            at = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
            # labels = c("low", "zero", "high"),
            title = "Dissimilarity",
            row_title = "I am a big column title",
            # font = 10,
            legend_height = unit(5, "cm"),
            grid_width = unit(0.5, "cm"),
            title_position = "leftcenter-rot"),
          na_col = "black",
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          row_split = factor(f, levels = unique(f)),
          column_split = factor(f, levels = unique(f)),
          width = unit(15, "cm"),
          height = unit(15, "cm"),
          row_gap = unit(1, "mm"),
          column_gap = unit(1, "mm"), border = FALSE)
}

f.corrplot <- function(m,
                       p,
                       type = "upper",
                       col = col.pallett){
  colnames(m) <- rownames(m) <- families %>% unique()
  corrplot(m,
           p.mat = p,
           method = "color",
           sig.level = 1e-10,
           type = type,
           is.corr = T,
           col = rep(col,2),
           # addCoef.col = col.pallett[1],
           diag = TRUE,
           cl.pos = "n",
           # tl.pos = "n",
           insig = "p-value")
  }

# correcting indeces ------------------------------------------------------

# df <- df.hhg.act  # create a copy of df
look <- snp.micro %>% select(Owner, ID, index)
# colnames(look) <- c("old.index", "new.index")
# look$old.index <- df %>% arrange(name.1) %>% pull(index.1) %>% unique() %>% as.numeric()
# look$new.index <- c(21:30,11:20,1:10,31:40,41:50)# %>% as.character()

## using lapply, loop over columns and match values to the look up table. store in "new".
# new <- df
# new$name.1 <- look$new.index[match(df$name.1, look$old.index)]
# new$name.2 <- look$new.index[match(df$name.2, look$old.index)]
# 
# 
# new.index <- look %>% arrange(old.index) %>% pull(new.index)
new.index <- snp.micro %>% arrange(fn) %>% pull(ID)
# new.index <- c(rev(new.index[1:10]),
#                rev(new.index[11:20]),
#                rev(new.index[21:30]),
#                rev(new.index[31:40]),
#                rev(new.index[41:50])
#                )
# lookup.replace <- function(t, look = look) lapply(t, function(x) look$new.index[match(x, look$old.index)])

df.hhg.act[,1:2] <- lapply(df.hhg.act[,1:2], function(x) look$index[match(x, look$Owner)])
df.hhg.act[,3:4] <- lapply(df.hhg.act[,3:4], function(x) look$ID[match(x, look$Owner)])

df.hhg.con[,1:2] <- lapply(df.hhg.con[,1:2], function(x) look$index[match(x, look$Owner)])
df.hhg.con[,3:4] <- lapply(df.hhg.con[,3:4], function(x) look$ID[match(x, look$Owner)])

# making Sim/Diff matrices ------------------------------------------------
# s.con <- s.con[rev(new.index),rev(new.index)]
# s.act <- s.act[rev(new.index),rev(new.index)]

dead.folks <- c("HD2",
                "HD3",
                "SC1",
                "SC3")

df.hhg.act[df.hhg.act$name.1 %in% dead.folks,5:12] <- NA
df.hhg.act[df.hhg.act$name.2 %in% dead.folks,5:12] <- NA
df.hhg.con[df.hhg.con$name.1 %in% dead.folks,5:12] <- NA
df.hhg.con[df.hhg.con$name.2 %in% dead.folks,5:12] <- NA

hhg.act <- df.hhg.act %>% make.hhg.mat()
hhg.con <- df.hhg.con %>% make.hhg.mat()

l.Dissimilarity <- list(`Anatomical HHG dissimilarity` = hhg.con,
                     `Anatomical NetSimile dissimilarity` = make.halved.mats(s.con, s.con, TRUE),
                     `Functional HHG dissimilarity` = hhg.act,
                     `Functional NetSimile dissimilarity` = make.halved.mats(s.act, s.act, TRUE)
                     ) #%>%
  # map(function(x) x[new.index,new.index])



l.Dissimilarities.halved <- list(`Anatomical` = make.hhg.s(l.Dissimilarity$`Anatomical HHG dissimilarity`,
                                                          l.Dissimilarity$`Anatomical NetSimile dissimilarity`),
                              `Functional` = make.hhg.s(l.Dissimilarity$`Functional HHG dissimilarity`,
                                                      l.Dissimilarity$`Functional NetSimile dissimilarity`)
                              )
l.Contrast <- l.Dissimilarity %>%
  map(~make.block.mean(make.halved.mats(.,.,TRUE), TRUE))

for(l in l.Contrast) rownames(l) <- colnames(l) <- unique(families)

Differentiation <- l.Contrast %>% 
  ldply(function(l){
    sm <- NULL
    for(i in 1:5){
      # sm <- 
      x <- l[i,]
      sm[i] <- (1-x[i])/(mean(1-x[-i]))
    }
    names(sm) <- families %>% unique()
    sm %>% return()
  }
  )

g.d <- Differentiation[c(2,4),] %>%
  mutate(.id = c("Anatomical", "Functional")) %>% 
  gather("Family", "Differentiation", -.id)


f <- families %>% unique() %>% rep(each = 2)
g.d$Family <- factor(f, levels = unique(f))


anatomical <- l.Dissimilarities.halved$Anatomical %>%
  sim.hm() %>% 
  draw() %>% 
  grid.grabExpr()

functional <- l.Dissimilarities.halved$Functional %>%
  sim.hm() %>% 
  draw(row_title = "Heatmap list") %>% 
  grid.grabExpr()

cowplot::plot_grid(anatomical, functional, nrow = 2)


colnames(g.d)[1] <- "Network"

f.corrplot(l.Contrast$`Anatomical NetSimile dissimilarity`, l.Contrast$`Anatomical HHG dissimilarity`)
f.corrplot(l.Contrast$`Functional NetSimile dissimilarity`, l.Contrast$`Functional HHG dissimilarity`)


# Plotting Dissimilarity/Contrast/Differentiation -------------------------

g.d %>% ggplot(aes(x = `Network`,
                   y = Differentiation,
                   fill = Family),
               title = "avg") +
  # ylim(-.42,.42) +
  geom_bar(position = position_dodge(width=0.85), size = 2, width = 0.80, stat="identity") +
  geom_hline(yintercept = 1,linetype="dashed")


# Plotting graphs ---------------------------------------------------------



## Change the index in this line
index <- 4 # should be 1:4, for NetSimile/HHG of Anat/Func
plot.name <- pull(Differentiation, .id)[index] %>% gsub(" .*","",.)
plot.name %>% paste("Graph -",.) %>% print()
differ <- Differentiation[index,] %>%
  select(-.id) %>% 
  as.numeric()

positive.differ <- rep("*", 5)
positive.differ[differ<1] <- ""
# positive.differ <- ""

contra <- l.Contrast[[index]]
colnames(contra) <- families %>%
  unique() %>% 
  paste0(positive.differ,.)

g <- contra %>% graph_from_adjacency_matrix(mode = "undirected",
                                            weighted = TRUE,
                                            diag = FALSE)

#Color scaling function
c_scale <- colorRamp(col.pallett)

#Applying the color scale to edge weights.
#rgb method is to convert colors to a character vector.
E(g)$width <- E(g)$weight
E(g)$color <- apply(c_scale(E(g)$weight), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
V(g)$color <- apply(c_scale(diag(contra)), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
# V(g)$size <- 100*abs(differ)

v.frame.color <- rep("green", 5)
v.frame.color[differ<1] <- "white"


radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}

la <- layout.circle(g)

par(mar=c(0,0,0,0)+1.7)
g %>% plot(edge.width = 40*E(g)$weight,
           main = gsub(" dissimilarity", "", plot.name),
           vertex.size = 40*abs(differ),
           vertex.label.dist = 3.5,
           # vertex.frame.color = v.frame.color,
           vertex.frame.width = 100,
           vertex.label.degree = radian.rescale(x=1:5, direction=-1, start=6),
           layout=la
           )



# compareing Contrast_Anatomical and Contrast_Functional ------------------
## Garbage!
# c.an <- l.Contrast$`Anatomical NetSimile dissimilarity` %>% 
#   make.halved.mats(.,.,T)
# c.fu <- l.Contrast$`Functional NetSimile dissimilarity`





