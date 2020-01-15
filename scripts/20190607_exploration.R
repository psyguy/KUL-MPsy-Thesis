m <- make_random_graph(size = num_nodes,
                              num.links = num_edges,
                              seed = seed)
a <- num_nodes %>% runif(-1, 1) %>% t()

c <-  m %>% my_clustceof()

l <- list(activities = a,
          mat.connectivity = list(m),
          coef.clustering = c)

p <- list(num_nodes = num_nodes,
          num_edges = num_edges)

toy_brain <- new("brain",
                 birthday = as.character(Sys.time()),
                 age = list(beat = 1, minute = 1),
                 starting_values = l,
                 parameters = p,
                 history = l,
                 now = l
                 )

toy_brain_1 <- toy_brain %>% update()














# graveyard ---------------------------------------------------------------


a <- matrix(
  c(
    0,1,0,0,0,0,0,
    1,0,1,0,1,0,0,
    0,1,0,1,1,0,0,
    0,0,1,0,0,1,0,
    0,1,1,0,0,1,0,
    0,0,0,1,1,0,1,
    0,0,0,0,0,1,0
  ),
  nrow = 7
)

tr(a)

a %*% a %*% a %>% tr()


a %>% graph_from_adjacency_matrix() %>% transitivity()

a %*% a %>% tr()

a



rm(list = ls())
source("./functions/functions_GongvLeeuwen2004.R")
seed <- 1
# number of nodes
num_nodes <- 30
# number of links
num_edges <- round(2 * log(num_nodes) * (num_nodes - 1))
# number of iterations
T_ <- 10
# tol <- 0.001


# start of the script -----------------------------------------------------

set.seed(seed)
conn <- make_random_graph(size = num_nodes,
                          num.links = num_edges,
                          seed = rnorm(1))
conn %>% my_clustceof()

conn %>% pimage()
o <- conn %>% seriate()
pimage(conn,o)



student(name="John", age=21, GPA=3.5)
