# rm(list=ls())
source("./functions/functions_extract.R")
load("./data-pc/last-snapshot-conn-activity_20191010_1926.RData")

#### NetSimile #################

library(gtools) 
library(sna)
library(igraph)
library(e1071)
library(HHG)

featExt <- function(A){
  # A <- adjvec %>% vec2mat()
  G <- graph_from_adjacency_matrix(A, "undirected")
  # A <- as.matrix(get.adjacency(G))
  v <- vcount(G)
  FM <- matrix(0, nrow=v, ncol=7) # Feature Matrix
  for(node in 1:v){
    #1) Degree: No of neighbours
    neighbours <- neighbors(G, node)
    FM[node,1] <- length(neighbours)
    
    #2) Clustering cofficient 
    FM[node,2] <- transitivity(G, type=c("local"), vids=node, isolates=c("zero"))
    
    #3) average number of node i's two hop away distance neighbours
    nsize <- neighborhood.size(G, order=1, nodes=neighbours)
    if(length(nsize) > 0){
      FM[node,3] <- mean(nsize)
    }else{
      FM[node,3] <- 0  # Replacing NaN with 0
    }
    
    # 4 average clustering cofficient of neighbours
    avgNC <- transitivity(G, type=c("local"), vids=neighbours, isolates=c("zero"))
    if(length(avgNC) > 0 ){
      FM[node,4] <- mean(avgNC)
    }else{
      FM[node,4] <- 0  # Replacing NaN with 0
    }
    
    # 5 No of edges in ego(i)
    FM[node,5] <- sapply(ego.extract(A, ego=node),sum)/2
    
    # 6 No of outgoing edges from ego(i)
    egoNetwork <- igraph::neighborhood(G, nodes=node, order=1)[[1]]
    allNeighbours <- unlist(igraph::neighborhood(G, order=1, nodes = egoNetwork))
    FM[node,6] <- length(allNeighbours[!allNeighbours %in% egoNetwork])
    
    # No of neighbours of ego(i)
    FM[node,7] <- sum(sapply(igraph::neighborhood(G, order=1, nodes=egoNetwork), length)) - length(egoNetwork)
  }
  FM
}

## Algo 3 Features Aggregation : Generating Signature Vectors 

featAggr <- function(FM){
  sig <- c()
  sig <- append(sig, apply(FM,2, median))    # Median of each feature
  sig <- append(sig, apply(FM,2, mean))      # Mean of each feature
  sig <- append(sig, apply(FM,2, sd, na.rm=T))        # Standard Deviation of each feature
  sig <- append(sig, apply(FM,2, skewness, na.rm=T))  # Skewness of each feature
  sig <- append(sig, apply(FM,2, kurtosis, na.rm=T))  # Kurtosis of each feature
  sig
}


calc_signature <- function(s) featAggr(featExt(s))


l.connectivities <- l.extracted %>% map("m")
l.activities <- l.extracted %>% map("a") %>% map(my_coherenceD)

#takes 130 + 392  seconds
signatures.connectivities <- l.connectivities %>% ldply(calc_signature)
signatures.activities <- l.activities %>% ldply(calc_signature)


# making feature distributions --------------------------------------------

system.time(features.connectivities <- l.connectivities %>% map(featExt))
system.time(features.activities <- l.activities %>% map(featExt))

# features.activities <- l.activities %>% ldply(featExt)

signatures.connectivities <- features.connectivities %>% ldply(featAggr)
signatures.activities <- features.activities%>% ldply(featAggr)

dist.signatures <- signatures.connectivities[-1] %>%
  dist(method = "canberra") %>% 
  as.matrix()

dist.signatures %>% pimage()

system.time(
  distances.connectivities <- features.connectivities %>% map(~as.matrix(dist(.x, diag = TRUE, upper = TRUE)))
            )

# writing nested for loop of HHG ------------------------------------------

i.and.j <- list()
k <- 1
for(i in 1:length(d)){
  for(j in (i):length(d)){
    i.and.j[[k]] <- c(i,j)
    k <- k + 1
  }
}

i <- i.and.j[[counter]][1]
j <- i.and.j[[counter]][2]

d <- distances.connectivities

hhg.results <- list(indices = c(i,j),
                    names = c(names(d[i]),names(d[j])),
                    hhg.values = hhg.test(dd[[i]],distances.connectivities[[j]])
                    )
save_vars("hhg.results", prefix = "HHG")
