# Notes -------------------------------------------------------------------

#
#   
#   cf. _trial. Added to that, the parameters now include vectors of
#   epsilon and a parameter, allowing partitioning.
#   
#
#
#
#
#


# init --------------------------------------------------------------------

rm(list = ls())

source("./functions/functions_trial.R")

vat_name <- NULL
vat_num_nodes <- 300
vat_num_edges <- 5200
vat_seed <- 2000
partitions <- list(eps = "0.3x.1, 0.4x.8, 0.5x.1",
                   a = "1.4x.4, 1.7x.2, 2x.4")

vat_eps <- make_paramvect(s = partitions$eps,
                          n = vat_num_nodes,
                          seed = vat_seed)

vat_a <- make_paramvect(s = partitions$a,
                        n = vat_num_nodes,
                        seed = vat_seed + 1)

parameters =  list(n_nodes = vat_num_nodes,
                   n_edges = vat_num_edges,
                   partitions = partitions,
                   eps = vat_eps,
                   a = vat_a,
                   seed = vat_seed,
                   lower_bound_starting = 0,
                   global_minmax = FALSE,
                   blind_swap = FALSE)

brain_case <- NULL
for(days in 1:10){
  brain_case <- trial_grow(parameters =  parameters,
                           n_rewires = 1000,
                           n_updates = 20,
                           freq_snapshot = 200,
                           name = vat_name,
                           brain_younger = brain_case,
                           quiet = FALSE)
  
  brain_case@history$coef.clustering %>% plot(main = paste(brain_case@name,
                                                           "at",
                                                           brain_case@age$rewires %/% 1000,
                                                           "days."))
  
}

save_vars(list.of.vars = "brain_case",
          prefix = paste0("brain_",
                          brain_case@name,
                          "_",
                          brain_case@parameters$n_edges,
                          "edges"))

