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

families <<- c("HC", "MC", "BL", "LC", "SC") %>% rep(each = 10)
# functions used for plotting heatmaps ------------------------------------

sim.hm <- function(m, title = "", hm_size = 15, col = col.pallett){
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
            # row_title = "I am a big column title",
            # font = 10,
            legend_height = unit(hm_size/3, "cm"),
            grid_width = unit(hm_size/30, "cm"),
            title_gp = gpar(fontfamily = "CMU Serif",
                            fontsize = 18),
            labels_gp = gpar(fontfamily = "CMU Serif",
                             fontsize = 18),
            title_position = "leftcenter-rot"),
          
          row_names_gp = gpar(fontfamily = "CMU Serif",
                              fontsize = 12),
          column_names_gp = gpar(fontfamily = "CMU Serif",
                                 fontsize = 12),
          
          row_title_gp = gpar(fontfamily = "CMU Serif",
                              fontsize = 20),
          column_title_gp = gpar(fontfamily = "CMU Serif",
                              fontsize = 20),
          na_col = "black",
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          row_split = factor(f, levels = unique(f)),
          column_split = factor(f, levels = unique(f)),
          width = unit(hm_size, "cm"),
          height = unit(hm_size, "cm"),
          row_gap = unit(hm_size/15, "mm"),
          column_gap = unit(hm_size/15, "mm"),
          border = FALSE)
}

plot_labels <- function(label = "Whole",
                        size = 5,
                        angle = 0,
                        family = "CMU Classical Serif"){
  ggplot() +
    annotate(geom = 'text',
             x=1, y=1,
             label = label,
             size = size,
             family = family,
             angle = angle) +
    theme_void()
}

f.corrplot <- function(m,
                       p,
                       type = "upper",
                       col = col.pallett){
  colnames(m) <- rownames(m) <- families %>% unique() # %>% paste(" ", ., " ")
  colnames(p) <- rownames(p) <- families %>% unique() # %>% paste(" ", ., "")
  par(family="CMU Serif", cex = 1.5)
  corrplot(m,
           p.mat = p,
           method = "color",
           sig.level = 1e-10,
           type = "upper",
           addgrid.col = "white",
           is.corr = T,
           col = rep(col,2),
           # addCoef.col = col.pallett[1],
           diag = TRUE,
           cl.pos = "n",
           tl.col = "black",
           tl.cex = 1.25,
           # tl.pos = "n",
           number.font = 90,
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

## First do not run the following 9 lines (that sets the dead networks to NA)
## and make figure 5, then come back and run them to calculate family measures.

dead.folks <- c("MC2",
                "MC3",
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


diag(l.Dissimilarities.halved$Anatomical) <- 0
hm_anatomical <- l.Dissimilarities.halved$Anatomical %>%
  sim.hm(hm_size = 20) %>% 
  draw() %>% 
  grid.grabExpr()

diag(l.Dissimilarities.halved$Functional) <- 0
hm_functional <- l.Dissimilarities.halved$Functional %>%
  sim.hm(hm_size = 20) %>% 
  draw() %>% 
  grid.grabExpr()

shelf(cowplot)

heatmaps <- plot_grid(hm_anatomical,
                      plot_labels(" "),
                      hm_functional,
                      nrow = 3,
                      rel_heights = c(1,0.05,1))

heatmaps_labels <- plot_grid(plot_labels("Anatomical", size = 15, angle = 90, family = "CMU Serif"),
                             plot_labels(" "),
                             plot_labels("Functional", size = 15, angle = 90, family = "CMU Serif"),
                             nrow = 3,
                             rel_heights = c(1,0.05,1))

fig5 <- plot_grid(heatmaps_labels,
                  heatmaps,
                  ncol = 2,
                  rel_widths = c(0.05,1))


save_plot("fig5.pdf",
          fig5,
          base_height = 1.6*11.69,
          base_width = 1.4*8.27)


colnames(g.d)[1] <- "Network"



f.corrplot(l.Contrast$`Anatomical NetSimile dissimilarity`, l.Contrast$`Anatomical HHG dissimilarity`)
f.corrplot(l.Contrast$`Functional NetSimile dissimilarity`, l.Contrast$`Functional HHG dissimilarity`)


col_fun = colorRamp2(seq(0,.69,0.1), col.pallett[1:7])
lgd = Legend(at = seq(0,0.7,0.1),
             # labels = c("low", "zero", "high"),
             title = "Contrast",
             col_fun = col_fun,
             # row_title = "I am a big column title",
             # font = 10,
             legend_height = unit(5, "cm"),
             grid_width = unit(1, "cm"),
             # direction = "horizontal", 
             # title_position = "topcenter"
             title_position = "leftcenter-rot",
             title_gp = gpar(fontfamily = "CMU Serif",
                             col = "Black", fontsize = 18),
             labels_gp = gpar(fontfamily = "CMU Serif",
                              col = "Black", fontsize = 18))
draw(lgd,
     x = unit(0.08, "npc"),
     y = unit(0.05, "npc"),
     just = c("left", "bottom"))



# Plotting Dissimilarity/Contrast/Differentiation -------------------------

fig7 <- g.d %>%
  ggplot(aes(x = `Network`,
             y = Differentiation,
             fill = Family),
             title = "avg") +
  # ylim(-.42,.42) +
  geom_bar(position = position_dodge(width=0.85),
           size = 2,
           width = 0.80,
           stat="identity") +
  theme_bw() +
  # labs(x = ifelse(fam.code == "SD", "Club size",element_blank()),
  #      y = ifelse(Partition_ == "whole", "Value",element_blank()),
  #      # title = paste("Normalized Rich-Club for",
  #      #               fam.code,
  #      #               paste0("(", Partition_, ")"))
  #      ) +
  theme(text = element_text(size = 20, family = "CMU Serif"),
        # axis.title.x = element_blank(),
        plot.title = element_text(size = 10,
                                  family = "CMU Serif",
                                  hjust = 0.5,
                                  vjust = -12),
        # plot.subtitle = element_text(size = 15,
        #                              family = font.cmu.serif.upright.italic,
        #                              hjust = 0.01,
        #                              vjust = -18),
        legend.position = "none"
        # axis.title.y = element_text("Differentiation",
        #                             size = 20, family = "CMU Serif")
  ) +
# scale_x_continuous(breaks = seq(0, max.clubsize, 20),
#                    limits = c(0, max.clubsize)) +
  scale_y_continuous(breaks = seq(0, 1.3, .25),
                     limits = c(0, 1.3)) +
  geom_text(aes(label = Family),
            size = 5,
            family = "CMU Serif",
            position = position_dodge(width=.85),
            vjust = -0.5) +
  # scale_x_discrete() +
  coord_fixed() +
  geom_hline(yintercept = 1,linetype="dashed") +
  scale_fill_viridis_d(direction = -1,
                       begin = 0,
                       end = 0.85) +
  xlab("Network") +
  ylab("Differentiation")

save_plot("fig7.pdf",
          fig7,
          base_height = 6,
          base_width = 8.27)


# Plotting graphs ---------------------------------------------------------



## Change the index in this line
index <- 2 # should be 2 or 4, for NetSimile of Anat/Func
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
# par(mar=c(0,0,0,0)+1.7)

# par(mfrow=c(1,1))

par(mar = c(0,0,0,10))
g %>% plot(edge.width = 30*E(g)$weight,
           # main = gsub(" dissimilarity", "", plot.name),
           # family = "CMU Serif",
           vertex.frame.color = NA, 
           vertex.label.family = "CMU Serif",
           vertex.label.cex = 2,
           vertex.size = 35*abs(differ),
           vertex.label.color = "black",
           # vertex.label.dist = 5.5,
           vertex.frame.color = NULL,
           vertex.frame.width = 100,
           # vertex.label.degree = radian.rescale(x=1:5, direction=-1, start=6),
           layout=la
           )

col_fun = colorRamp2(seq(0,.79,0.1), col.pallett[1:8])
lgd = Legend(at = seq(0,0.8,0.2),
             # labels = c("low", "zero", "high"),
             title = "Contrast",
             col_fun = col_fun,
             # row_title = "I am a big column title",
             # font = 10,
             legend_height = unit(5, "cm"),
             grid_width = unit(1, "cm"),
             # direction = "horizontal", 
             # title_position = "topcenter"
             title_position = "leftcenter-rot",
             title_gp = gpar(fontfamily = "CMU Serif",
                             col = "Black", fontsize = 18),
             labels_gp = gpar(fontfamily = "CMU Serif",
                              col = "Black", fontsize = 18))
draw(lgd,
     x = unit(0.8, "npc"),
     y = unit(0.3, "npc"),
     just = c("left", "bottom"))


## save the resulting plot at 6 x 9 inches



# x <- recordPlot()
# 
# title(gsub(" dissimilarity", "", plot.name),
#       family = "CMU Serif",
#       cex.main=3,
#       col.main="green")
# 
# plot_grid(x,x, nrow = 2)


# compareing Contrast_Anatomical and Contrast_Functional ------------------
## Garbage!
# c.an <- l.Contrast$`Anatomical NetSimile dissimilarity` %>% 
#   make.halved.mats(.,.,T)
# c.fu <- l.Contrast$`Functional NetSimile dissimilarity`





