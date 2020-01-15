source("./scripts/20191202_family-comparison-heatmaps.R")

library(plyr)
library(viridis)

hm <- function(m, title = ""){
  f <- families
  if(nrow(m)==5) f <- families %>% unique()
  m <- m %>% as.matrix() %>% unname()
  Heatmap(m,
          col = inferno(10,1,0,1,-1),
          #colorRamp2(c(0, 1),c("white", "azure4")),
          # col = (rainbow(100, start = 0, end = 1)),
          column_title = title,
          name = "Dissimilarity",
          na_col = "white",
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.1f", small_mat[i, j]), x, y, gp = gpar(fontsize = 10))
          },
          # column_names_rot = -45, column_names_side = "top", right_annotation = ha,
          row_split = factor(f, levels = unique(f)),
          column_split = factor(f, levels = unique(f)),
          row_gap = unit(0.5, "mm"), column_gap = unit(0.5, "mm"), border = TRUE)

}


l.family.halved[[1]] %>% hm()

l <- l.avg[[2]]
rownames(l) <- colnames(l) <- unique(families)

for(l in l.avg) rownames(l) <- colnames(l) <- unique(families)

m <- l
f.corrplot <- function(m, type = "lower") corrplot(m,
                                                   method = "shade",
                                                   type = type,
                                                   # title = "\n Connectivity HHG test \n",
                                                   is.corr = T,
                                                   col = rep(inferno(10,1,0,1,-1),2),
                                                   # addgrid.col = NA,
                                                   # addCoefasPercent = T,
                                                   addCoef.col = "black",
                                                   diag = TRUE,
                                                   # tl.pos = "ld",
                                                   cl.pos = "n",
                                                   tl.pos = "n",
                                                   # tl.col = "white",
                                                   number.cex = .6) #%>% png()

l.avg$`Connectivity HHG test` %>% f.corrplot()


d <- l.family.sim$`Connectivity HHG test`
d[1:10,11:50] %>% sum

f <- families %>% unique()
sm <- NULL


for(i in 1:5){
  x <- l.avg[[1]][i,]
  sm[i] <- (x[i])/(mean(x[-i]))
}
names(sm) <- f

group.differentiation <- l.avg %>% 
  ldply(function(l){
    for(i in 1:5){
      x <- l[i,]
      # sm[i] <- (10-x[i])/(sum(50-x[-i]))
      sm[i] <- (1-x[i])/(mean(1-x[-i]))
    }
    names(sm) <- f
    sm %>% return()
    }
    )

g.diff <- function(m, dis){
  m[is.na(m)] <- 0
  n <- max(m)-m
  # n <- m
  x <- NULL
  for(i in 1:5){
    n.tmp <- n[(10*i-9):(i*10),] %>% sum(na.rm = TRUE)
    n.nume <- n[(10*i-9):(i*10),(10*i-9):(i*10)] %>% sum(na.rm = TRUE)
    x[i] <- n.nume/(n.tmp-n.nume)
  }
  names(x) <- families %>% unique()
  x %>% return()
}

# group.differentiation <- l.family.sim %>% 
#   ldply(g.diff)

g.d <- group.differentiation %>% gather("Family", "Differentiation", -.id)

g.d$Family <- factor(g.d$Family, levels = unique(g.d$Family))

colnames(g.d)[1] <- "Resemblance measure"
g.d %>% ggplot(aes(x = `Resemblance measure`, y = Differentiation, fill = Family), title = "avg") +
  geom_bar(position="dodge", stat="identity")


