# rm(list = ls())
source("./functions/functions_extract.R")

load("./data-pc/snp_only-1e+6_20190723_1404.RData")

m <- snp$adj.mat.vect[[1]] %>% vec2mat()


# SVD and eigen stuff -----------------------------------------------------


o.eigendecom <- miRNAss::eigenDecomposition(m,1)#$U
# o.eigen <- eigen(m)
s <- m %>% seriate()
ser <- m[s[[1]], s[[2]]]
o.svd <-  ser %>% svd()


o.eigendecom <- miRNAss::eigenDecomposition(ser,1)#$U


e.v <- matrix(0,300,300)
e.v[1,] <- o.eigendecom$U
mid <- matrix(0,300,300)
mid[1,1] <- o.eigendecom$D


mid <- diag(o.svd$d)
red <- 300:300
left <- o.svd$u
left[red,] <- 0
left <- left %>% t()

(x <- left %*% mid %*% t(left)) %>%
  pimage(key = FALSE, main = "kill me")



(x <- e.v %*% mid %*% t(e.v)) %>%
  pimage(key = FALSE, main = "kill me")

d.red[red, red] <- 0



(m.rec <- o.svd$v %*% d.red) %>% 
  pimage(key = FALSE, main = "seriated")


o.svd$u %>% pimage(key = FALSE)
o.svd$v %>% pimage(key = FALSE)
m.rec %>% pimage(key = FALSE, main = "reduced")
ser %>% pimage(key = FALSE, main = "seriated")

(m.rec <- o.svd$v - (o.svd$u)) %>% 
  pimage()#key = FALSE, main = "seriated")

sum(m-m.rec)


svdmaker <- function(v){
  m <- v %>% vec2mat()
  s <- m %>% seriate()
  eig <- m[s[[1]], s[[2]]] %>% svd()
  eig$d %>% return()
}

eig$v %>% pimage


svd.e <- eig$d

p <-  m[s[[1]], s[[2]]] %>% princomp()

p$loadings




snp.whole <- snp %>% filter(Partition=="whole")

e <- snp.whole$adj.mat.vect %>% map(svdmaker)

num.persons <- 50
eigensim <- matrix(0, num.persons, nrow = num.persons)
for(i in 1:num.persons){
  for(j in 1:num.persons){
    eigensim[i,j] <- e[[i]] %*% e[[j]]
  }
}

e.centered <- eigensim-(max(eigensim)+min(eigensim))/2
es <- e.centered %>% seriate()

pimage(e.centered, key=F)
pimage(e.centered, es, key=F)

ee <- e %>% simplify2array %>% t()

k <- kmeans(ee, 5)
k.seriated <- 

k$cluster %>% hist()

cl <- k$cluster %>% cbind(snp.whole$Verbal.Description %>% as.character())

ee

h <- ee %>%
  dist(method = "euclidean") %>% 
  hclust()

h %>% plot()






# RV and MatrixCorrelation ------------------------------------------------


library("MatrixCorrelation")

m1 <- snp$adj.mat.vect[[1]] %>% vec2mat()
m2 <- snp$adj.mat.vect[[121]] %>% vec2mat()

MatrixCorrelation::allCorrelations(m1,m1)
nc <- 300
MatrixCorrelation::allCorrelations(m1,m2,nc,nc)

red <- 20:300
X1  <- m1# scale( matrix( rnorm(100*300), 100,300), scale = FALSE)
usv <- svd(X1)
X2  <- usv$u[,-red] %*% diag(usv$d[-red]) %*% t(usv$v[,-red])
PSI(X1,X2)
nc <- min(red)-1
MatrixCorrelation::allCorrelations(X1,X2)#,nc,nc)
