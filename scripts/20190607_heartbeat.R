# Notes -------------------------------------------------------------------

#
#   
#   at each heartbeat nodes are updated (logistic functions)
#   and every minute the rewiring takes place
#
#
#
#
#
#


# init --------------------------------------------------------------------

rm(list = ls())
source("./functions/functions_heartbeat.R")
seed <- 1

# number of nodes
num_nodes <- 300
# number of links
num_edges <- round(2 * log(num_nodes) * (num_nodes - 1))

# number of heartbeats per minute
num_hr <- 10
# number of rewirings
num_minutes <- 500


# start of the script -----------------------------------------------------

set.seed(seed)

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


aging_brain <- toy_brain


time_start <- Sys.time()

for(i in 1:100) {
  aging_brain <- aging_brain %>% update()
  if (!(i %% 10))
    aging_brain <- aging_brain %>% rewire()
}

time_taken <- Sys.time - time_start

save_vars(ls())

# GvL.finished$connectivity.matrix %>% corrplot(method = "square")

# GvL.finished$x.tot[, ]
#
# (GvL.finished$x.tot[2, ] - GvL.finished$x.tot[3, ]) %>% sum()


# for(k in c(1:11)){

# densityplot(GvL.finished$x.tot[50,])
#
# GvL.finished$x.tot %>% corrplot()




