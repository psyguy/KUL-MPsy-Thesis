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

alpha <- c(0.5, 1, 1, 1, 2, 5)
beta <- c(0.5, 1, 2, 5, 1, 1)
l_ <- data.frame(alpha,beta)

vat_name <- NULL
vat_num_nodes <- 300
vat_num_edges <- 5200
round <- 1
vat_seed <- 1639

row_eps <- 3
row_a <- 4

params.eps_a <- l_[row_eps,] %>% cbind(l_[row_a,]) %>%
  format(nsmall=1)
colnames(params.eps_a) <- c("eps", "eps", "a", "a") %>%
  paste(colnames(params.eps_a), sep = ".")


vat_eps <- make_paramdist(alpha_beta = params.eps_a[1:2],
                          range_param = c(0.3,0.5),
                          n = vat_num_nodes,
                          seed = vat_seed)

vat_a <- make_paramdist(alpha_beta = params.eps_a[3:4],
                        range_param = c(1.4,2),
                        n = vat_num_nodes,
                        seed = vat_seed + 1)

parameters =  list(params.eps_a = params.eps_a,
                   round = round,
                   n_nodes = vat_num_nodes,
                   n_edges = vat_num_edges,
                   seed = vat_seed,
                   round = round,
                   lower_bound_starting = 0,
                   brain.code <- "",
                   eps = vat_eps,
                   a = vat_a)

# brain_case <- NULL
for(days in 1:1){
  brain_case <- trial_grow(parameters =  parameters,
                           n_rewires = 1000,
                           n_updates = 20,
                           freq_snapshot = 200,
                           name = vat_name,
                           brain_younger = brain_case,
                           quiet = FALSE)
  
  # brain_case@history$coefficients %>%
  #   plot(main = paste(brain_case@name,
  #                     "at",
  #                     brain_case@age$rewires %/% 1000,
  #                     "days."))
  
}
# 
# save_vars(list.of.vars = "brain_case",
#           prefix = paste(brain_case@parameters$brain.code, brain_case@name, sep = "_"))


