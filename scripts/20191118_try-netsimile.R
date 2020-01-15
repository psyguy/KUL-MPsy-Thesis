# rm(list=ls())
source("./functions/functions_extract.R")


#### NetSimile #################

library(gtools) 
library(sna)
library(igraph)
library(e1071)

## Generating igraph instances from individual graph files
getGraph <- function(x){ 
  x <- read.table(x)         # Reading graph from the directory
  if(x[1,2] <= 0){           # Checking for empty edges
    return()
  }
  gv <- x[1,1]   						 # Fetching the vertex count	
  ge <- x[2:nrow(x),] + 1    # Adding 1 to move vertex id, because of 0 vertex-id
  g <- graph.empty() + vertices(1:gv)
  g <- add.edges(g, t(ge))
  g <- as.undirected(g)      # making the graph undirected
  g
}

## Algo 2 Extracting features from the graph
featExt <- function(adjvec){
  A <- adjvec %>% vec2mat()
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
  sig <- append(sig, apply(FM,2, sd))        # Standard Deviation of each feature
  sig <- append(sig, apply(FM,2, skewness, na.rm=T))  # Skewness of each feature
  sig <- append(sig, apply(FM,2, kurtosis, na.rm=T))  # Kurtosis of each feature
  sig
}

##### Algo 1: NetSimile, computing similarity scores between pair of graphs

adjvect.1 <- snp$adj.mat.vect[1]
adjvect.2 <- snp$adj.mat.vect[5]
adjvect.3 <- snp$adj.mat.vect[101]

s <- adjvect.1 %>% featExt()
s %>% featAggr()
s[,7] %>% densityplot()

sigcalc <- function(s) featAggr(featExt(s))

adj.list <- list(adjvect.1, adjvect.2, adjvect.3)

sigvectors <- ldply(adj.list, sigcalc)

d <- sigvectors %>% dist(method = "canberra")

d %>% hclust() %>% plot()


system.time(
sigvectors <- snp %>% 
  filter(Partition == "whole") %>% 
  filter(Rewiring == 1e6) %>% 
  # filter(Verbal.Description == "Hyper-coupled minority") %>% 
  pull(adj.mat.vect) %>%
  ldply(sigcalc)
)


for(i in 1:5){
  
}
d3 <- sigvectors[c(21:30),] %>%
  dist(method = "canberra") %>%
  as.vector() #%>% 
  densityplot()

di %>% pimage

di <- sigvectors %>%
  dist(method = "canberra") # %>% hclust() %>% plot()
c.d <- cluster::diana(di)
pltree(c.d, cex = 0.6, hang = -1,
       main = "RV")
rect.hclust(c.d, k = 5, border = 2:5)


di %>% hclust(method = "ward.D") %>% plot()


wilcox.test(as.numeric(sigvectors[1,]),as.numeric(sigvectors[4,]))




netsimile <- function(){
  sigVectors <- list()
  for(i in 1:3){
    file <- paste(directory, fileNames[i], sep="/")
    if(file.exists(file)){
      graph <- getGraph(file)             # Generating igraph instance
      if(is.null(graph)){  next } 
      features <- featExt(graph)          # Extracting feature matrix (Algo1)
      sigVectors <- append(sigVectors, list(featAggr(features)))  # Aggregation (Algo2)
      cat("  Signature Vector computed for file :", file, "\n")  
    }
  }
  cat("  Signature Vectors for the graphs has been calculated \n")
  
  # Computing Distance measure as scores using canberra distance 
  # between adjacent 
  scores <- c()
  for(j in 1:length(sigVectors)){
    if((j+1) <= length(sigVectors)){
      scores <- append(scores, dist(rbind(sigVectors[[j]], sigVectors[[j+1]]), method="canberra"))        
    }
  }
  
}  
