source("./functions/functions_extract.R")

load("./data-pc/snp_only-1e+6_20190723_1404.RData")

m <- snp$adj.mat.vect[[1]] %>% vec2mat()

# m200 <- snp$adj.mat.vect[[200]] %>% vec2mat()
# m100 <- snp$adj.mat.vect[[100]] %>% vec2mat()
# 
# e200 <- eigen(m200)$values
# e100 <- eigen(m100)$values
# 
# e200 %*% e100

# matsim <- function(v1, v2){
#   
#   m1 <- v1 %>% vec2mat()
#   m2 <- v2 %>% vec2mat()
#   
#   s1 <- m1 %>% seriate()
#   s2 <- m2 %>% seriate()
#   
#   eig1 <- m1[s1[[1]],s1[[2]]] %>% eigen(symmetric = TRUE)
#   eig2 <- m2[s2[[1]],s2[[2]]] %>% eigen(symmetric = TRUE)
#   
#   sim <- eig1$values %*% eig2$values
#   sim %>% return()
# }
# 
# snpmatsim <- function(r1, r2, snp){
#   matsim(snp$adj.mat.vect[[r1]],snp$adj.mat.vect[[r2]])
# }
# 
# eigenmaker <- function(v){
#   
#   m <- v %>% vec2mat()
# 
#   s <- m %>% seriate()
# 
#   eig <- m[s[[1]], s[[2]]] %>% eigen(symmetric = TRUE)
#   
#   eig$values %>% return()
# }


svdmaker <- function(v){
  m <- v %>% vec2mat()
  s <- m %>% seriate()
  eig <- m[s[[1]], s[[2]]] %>% svd()
  eig$d %>% return()
}

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

# 
# #######################
# #Spectral Clustering
# ######################
# #Gaussian mixtures generation
# library(MASS)#load MASS library
# Sigma = diag(2)/10#Diagonal covariance matrix
# data1 = mvrnorm(100,c(0,0),Sigma)#group 1 with mean (0,0)
# data2 = mvrnorm(100,c(2,0),Sigma)#group 2 with mean (2,0)
# data3 = mvrnorm(100,c(0,4),Sigma)#group 3 with mean (0,4)
# #Concatenate 3 gaussian clusters
# data = rbind(data1,data2,data3)
# #Linear covariance matrix computation
# scal_matrix = data%*%t(data)
# Kmatrix = as.matrix(dist(data,method='euclidean'))
# #Threshold and construction of the similarity-based adjacency matrix
# epsilon = 0.5
# A = (Kmatrix*epsilon)*1
# #Igraph construction based on the Adjacency matrix
# g = graph_from_adjacency_matrix(A,mode = 'undirected')
# #Graph Laplacian computation
# degrees = apply(A, 1,sum)
# L = diag(degrees)/300-A
# #Spectral clustering : 5 lines of code
# eigvalues = eigen(L)$values#SVD decomposition
# eigV =eigen(L)$vectors#Keep eigenvectors
# eigV = eigV [,300:1]#re-ordering of the eigenvectors
# U = eigV[,1:10]#Projection
# xx <- kmeans(U,4)#Clustering
