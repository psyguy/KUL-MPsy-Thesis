# Notes -------------------------------------------------------------------

#
#   
#   at each heartbeat nodes are updated (logistic functions)
#   and every minute the rewiring takes place
#   now brains have names, and eps = .3 (instead of .8!)
#
#
#
#
#


# init --------------------------------------------------------------------

rm(list = ls())

source("./functions/functions_trial.R")
load("./data/outcomes_20190612_1952.RData")
load("./data/brain_aging_20190612_1955.RData")

outcomes_decent <- outcomes %>%
  filter(mean_variance>.1) %>%
  arrange(coef.clustering_normalized) %>% 
  tail(10)

case_name <- "Emma Segers"
case_params <- outcomes_decent %>% filter(name == case_name)
# the best brain is Emma Segers. Simulating her for a longer time

brain_case <- trial_grow(parameters =  list(n_nodes = case_params$num_nodes,
                                          n_edges = case_params$num_edges,
                                          eps = case_params$eps,
                                          seed = case_params$seed),
n_rewires = 100,
n_updates = 20,
freq_snapshot = 200,
name = case_params$name,
save_brain = FALSE,
quiet = FALSE)

save_vars() #"brain_case", prefix = "brain_case")
