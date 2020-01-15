# Notes -------------------------------------------------------------------

# 
# It turns out that the #edges << N^2. And thats the reason why the 
# nothing changed over time. Will correct in the next version. 
# 
# 
# 
# 
# 
# 
# 
# 


# init --------------------------------------------------------------------

rm(list = ls())
source("./functions/functions_GongvLeeuwen2004.R")
seed <- 1
# number of links
L_c <- 520
# number of nodes
N <- 40#3000 #3000
# number of iterations
T_ <- 600 #600
# tol <- 0.001


# start of the script -----------------------------------------------------

set.seed(seed)

conn <- make_random_graph(size = N,
                          num.links = L_c,
                          seed = seed)

x.init <- N %>% runif(-1, 1) %>% t()

g <- igraph::graph_from_adjacency_matrix(conn, mode = "undirected")
ClCoef <- igraph::transitivity(g)

GvL.start <- list(
  x.tot = x.init,
  clustering.coefficient = ClCoef,
  connectivity.matrix = conn
)

GvL.finished <- GvL.start

for (i in 1:10) {
  print(i)
  GvL.finished <- GvL.finished %>%
    GongvLeeuwen2004_adaptive_rewire()
  g <- GvL.finished$connectivity.matrix %>%
    igraph::graph_from_adjacency_matrix(mode = "undirected")
  (ClCoef <- igraph::transitivity(g)) %>% print()
}

GvL.finished$x.tot[, ]

(GvL.finished$x.tot[2, ] - GvL.finished$x.tot[3, ]) %>% sum()


# for(k in c(1:11)){

densityplot(GvL.finished$x.tot[1,])

