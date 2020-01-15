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

# case_name <- "Emma Segers"
# case_params <- outcomes_decent %>% filter(name == case_name)

# making Hellrigel's brain
case_params <- outcomes_decent[1,]
case_params[1,] <- NA
case_params$name <- NULL#"Stefan Hellrigel 2Xedges"
case_params$eps <- 0.8
case_params$num_nodes <- 300
case_params$num_edges <- 5200
case_params$seed <- 1159

parameters =  list(n_nodes = case_params$num_nodes,
                   n_edges = case_params$num_edges,
                   eps = case_params$eps,
                   seed = case_params$seed,
                   lower_bound_starting = 0,
                   global_minmax = FALSE,
                   blind_swap = FALSE)

# brain_case <- NULL
for(days in 1:10){
  brain_case <- trial_grow(parameters =  parameters,
                           n_rewires = 1000,
                           n_updates = 20,
                           freq_snapshot = 200,
                           name = case_params$name,
                           brain_younger = brain_case,
                           quiet = FALSE)
  
  brain_case@history$coef.clustering %>% plot(main = paste(brain_case@name,
                                                           "at",
                                                           brain_case@age$rewires %/% 1000,
                                                           "days."))
  
  
  tmp_seq <- seq(1, brain_case@age$rewires, 10)
  
  sbs_h <- (brain_case@history$activities[tmp_seq,])
  
  sbs_h[,12] %>% plot()
  
  sbs_h %>% gplots::heatmap.2(dendrogram = 'none',
                              Rowv = FALSE,
                              Colv = FALSE,
                              margins = c(1, 1),
                              col = colorRampPalette(c("white","yellow","orange","red"))(n = 299),#brewer.pal(name = "RdBu"),
                              # key = FALSE,
                              # density.info = "none",
                              trace = 'none',
                              xlab = "nodes",
                              ylab = "rewirings",
                              main = paste(brain_case@name,
                                           "\n with",
                                           brain_case@parameters$n_edges,
                                           "edges \n after",
                                           brain_case@age$rewires-1,
                                           "rewirings")
  )
  
  
}
# save_vars()
save_vars(list.of.vars = "brain_case",
          prefix = paste0("brain_",
                          brain_case@name,
                          "_",
                          brain_case@parameters$n_edges,
                          "edges"))



# looking at the activations ----------------------------------------------
brain_case@name
brain_case@age$rewires
brain_case@parameters

Sys.time()
