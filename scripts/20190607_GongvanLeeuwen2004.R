# Notes -------------------------------------------------------------------

#
#   relation between #nodes and #edges:
#   similar to Jarman et al. 2017.
#
#   Ilias:N = 100, T_ = 14000
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
# number of nodes
num_nodes <- 30
# number of links
num_edges <- round(2 * log(num_nodes) * (num_nodes - 1))
# number of iterations
T_ <- 20000 #600
# tol <- 0.001


# start of the script -----------------------------------------------------

set.seed(seed)

conn <- make_random_graph(size = num_nodes,
                          num.links = num_edges,
                          seed = seed)

x.init <- num_nodes %>% runif(-1, 1) %>% t()

g <- igraph::graph_from_adjacency_matrix(conn, mode = "undirected")
ClCoef <- igraph::transitivity(g)

GvL.start <- list(
  x.tot = x.init,
  clustering.coefficient = ClCoef,
  connectivity.matrix = conn
)

GvL.finished <- GvL.start

time_start <- Sys.time()
# T_ <- 4000 #since more is redundant
for (i in 1:T_) {
  print(i)
  GvL.finished <- GvL.finished %>%
    GongvLeeuwen2004_adaptive_rewire()
  g <- GvL.finished$connectivity.matrix %>%
    igraph::graph_from_adjacency_matrix(mode = "undirected")
  (ClCoef <- igraph::transitivity(g)) %>% print()
  if (!(i %% (T_ / 40))) {
    title <- paste0("t = ", i, ", Clustering Coefficient = ", ClCoef)
    d <- GvL.finished$connectivity.matrix
    o <- seriate(d)
    pimage(d, o, main = title)
    # GvL.finished$connectivity.matrix %>% corrplot(method = "square")
  }
  
}


time_taken <- Sys.time - time_start
GvL.finished$clustering.coefficient %>% plot()

save_vars(ls())

# GvL.finished$connectivity.matrix %>% corrplot(method = "square")

# GvL.finished$x.tot[, ]
#
# (GvL.finished$x.tot[2, ] - GvL.finished$x.tot[3, ]) %>% sum()


# for(k in c(1:11)){

# densityplot(GvL.finished$x.tot[50,])
#
# GvL.finished$x.tot %>% corrplot()
