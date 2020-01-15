# rm(list = ls())
source("./functions/functions_extract.R")

install.packages("usedist")
library(factoextra)
library(cluster)
load("./data-pc/snp_only-1e+6_20190723_1404.RData")

m <- snp$adj.mat.vect[[1]] %>% vec2mat()


# SVD and eigen stuff -----------------------------------------------------


# RV and MatrixCorrelation ------------------------------------------------

my_seriate <- function(v){
  m <- ifelse(is.matrix(v), v, vec2mat(v))
  s <- m %>% seriate()
  m[s[[1]], s[[2]]] %>% return()
}


m1 <- snp$adj.mat.vect[[1]] %>% my_seriate()
m2 <- snp$adj.mat.vect[[121]] %>% my_seriate()

MatrixCorrelation::allCorrelations(m1,m1)
nc <- 300
x <- MatrixCorrelation::allCorrelations(m1,m2, methods = "RVadj")#,nc,nc)

red <- 20:300
X1  <- m1# scale( matrix( rnorm(100*300), 100,300), scale = FALSE)
usv <- svd(X1)
X2  <- usv$u[,-red] %*% diag(usv$d[-red]) %*% t(usv$v[,-red])
PSI(X1,X2)
nc <- min(red)-1
MatrixCorrelation::allCorrelations(X1,X2)#,nc,nc)

s <- snp %>%
  filter(Partition == "whole") %>% 
  pull(adj.mat.vect)
  
usedist::dist_make(s, MatrixCorrelation::RVadj, "RVadj")

num.persons <- 50
sim <- matrix(0, num.persons, nrow = num.persons)
system.time(
for(i in 1:num.persons){
  m1 <- s[[i]] %>% my_seriate()
  print(i)
  for(j in i:(num.persons)){
    m2 <- s[[j]] %>% my_seriate()
    sim[j,i] <- 1 - MatrixCorrelation::PSI(m1,m2)
    print(j)
  }
}
)

RVdist <- (1 - RVsim) %>% as.dist()
dist.PSI <- (1 - sim) %>% as.dist()

dist.PSI %>% hclust() %>% plot()

c <- dist.PSI %>% hclust(method = "ward.D")


di <- RVdist
c.d <- cluster::diana(di)
pltree(c.d, cex = 0.6, hang = -1,
       main = "RV")
rect.hclust(c.d, k = 5, border = 2:5)


fviz_dist(dist.PSI, gradient = list(low = "white", high = "black"))
