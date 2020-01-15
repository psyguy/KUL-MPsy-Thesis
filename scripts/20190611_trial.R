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
seed <- -99

set.seed(seed)
eps <- seq(.25, .45, by = 0.05)
num_nodes <- seq(300, 1000, by = 100)
mult_edges <- seq(1, 3, 0.5)

outcomes <- data.frame(matrix(nrow = 0, ncol = 9))
colnames(outcomes) <- c("name", "age_beat", "age_minute",
                      "num_nodes", "num_edges", "eps",
                      "seed", "mean_variance", "coef.clustering_normalized")



# looping over different parameterizations --------------------------------

time_start <- Sys.time()
for(i_eps in eps){
  for (i_num_nodes in num_nodes) {
    for(i_mult_edges in mult_edges){
      num_edges <- round(i_mult_edges * 2 * log(i_num_nodes) * (i_num_nodes - 1))
      if(i_num_nodes*(i_num_nodes-1)/2 < num_edges/2) next()
      paste("Now doing for epsilon =", i_eps, "and", i_num_nodes, "nodes &", num_edges, "edges") %>%
        print()
      toy_brain <- trial_vat(parameters =  list(num_nodes = i_num_nodes,
                                                num_edges = num_edges,
                                                eps = i_eps,
                                                freq_snapshot = 20,
                                                seed = 1000 * round(rnorm(1),3) %>% abs()
                                                  ),
                             num_minutes = 20,
                             num_hr = 20,
                             quiet = T)
      outcomes <- outcomes %>% rbind(trial_summary(toy_brain))
    }
  }
}
(time_taken <- Sys.time() - time_start)



outcomes$name <- outcomes$name %>% as.character()
indx <- sapply(outcomes, is.factor)
outcomes[indx] <- lapply(outcomes[indx], function(x) as.numeric(as.character(x)))
rm(indx)


save_vars(prefix = "trials")


# deciding which parameters to keep ---------------------------------------

outcomes_decent <- outcomes %>%
  filter(mean_variance>.1) %>%
  arrange(coef.clustering_normalized) %>% 
  tail(10)

case_name <- "Emma Segers"
case_params <- outcomes_decent %>% filter(name == case_name)
# the best brain is Emma Segers. Simulating her for a longer time

toy_brain <- trial_vat(parameters =  list(num_nodes = case_params$num_nodes,
                                          num_edges = case_params$num_edges,
                                          eps = case_params$eps,
                                          freq_snapshot = 50,
                                          seed = case_params$seed
                                      ),
                                      num_minutes = 2000,
                                      num_hr = 20,
                                      name = case_params$name,
                                      save_brain = TRUE,
                                      quiet = FALSE)
