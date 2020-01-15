# Intro -------------------------------------------------------------------

## This script keeps all functions needed for GongvLeeuwen2004.
## It's standalone and loads functions_my automatically

# rm(list = ls())


# loading the basic functions and packages --------------------------------

if(!exists("path.to.functions_my")) path.to.functions_my <- "functions"
source(paste0(path.to.functions_my,"/functions_my.R"))
source(paste0(path.to.functions_my,"/functions_heartbeat.R"))
source(paste0(path.to.functions_my,"/functions_netmeas.R"))



# a unified function to grow (update+rewire) ------------------------------

trial_grow <- function(parameters, # = NULL,
                       # parameters =  list(params.eps_a = c(0.5,0.5,1,2),
                       #                    round = 0,
                       #                    n_nodes = 100,
                       #                    n_edges = 0,
                       #                    seed = -99,
                       #                    round = 0,
                       #                    lower_bound_starting = 0,
                       #                    brain.code <- NULL,
                       #                    eps = NULL,
                       #                    a = NULL),
                       n_updates = 1, # number of heartupdates per rewireupdates of notes
                       n_rewires = 1, # number of rewirings of the network
                       freq_snapshot = 200, # frequency of saving connectivity matrices in the brain
                       brain_younger = NULL,
                       save_brain = FALSE,
                       name = NULL,
                       quiet = TRUE){
  
  
  brain_growing <- brain_younger

  
# making a brain if there is no younger brain to continue growing ---------

  if(is.null(brain_growing)){
    
    # if(is.null(parameters)) stop("Gimme the parameter list!")
    if(!parameters$n_edges){
      parameters$n_edges <- round(1.5 * 2 * log(parameters$n_nodes) * (parameters$n_nodes - 1))
    }
    if(parameters$seed == 0 | parameters$seed == -99){
      if(parameters$seed == -99) set.seed(-99)
      parameters$seed <- 1000 * round(rnorm(1),3) %>% abs()
    }
    

    
    # making the name, seed, and brain code
    bc_ <- make_brain.code(parameters)
    parameters$brain.code <- bc_$braincode
    parameters$seed <- bc_$seed
    name <- bc_$name
    
    m <- make_random_graph(size = parameters$n_nodes,
                           num.links = parameters$n_edges,
                           seed = parameters$seed)
    act <- parameters$n_nodes %>%
      runif(parameters$lower_bound_starting, 1) %>%
      t()
    
    # setting eps and a parameter vectors

    l <- list(
      activities = act,
      mat.connectivity = list(m),
      coefficients = NULL
    )
    
    
    brain_growing <- new(
      "brain",
      name = name,
      birthday = as.character(Sys.time()),
      age = list(updates = 1, rewires = 1),
      parameters = parameters,
      initial = l,
      history = l,
      now = l
    )
    
    brain_growing@initial$coefficients <-
      brain_growing@now$coefficients <- NULL
        # netmeas_coefs(b = brain_growing, t_ = 1, name = name)
    
  }
  
# writing the update and rewire again together here -----------------------


  # extracting values of now
  m.l <- brain_growing@now$mat.connectivity -> now.m
  if(is.list(m.l)) now.m <- m.l[[length(m.l)]]
  now.a <- brain_growing@now$activities
  # now.c <- brain_growing@now$coef.clustering
  eps <- brain_growing@parameters$eps
  a <- brain_growing@parameters$a

  now.updates <- brain_growing@age$updates
  now.rewires <- brain_growing@age$rewires
  
    
  # extracting historical values
  his.a <- brain_growing@history$activities
  his.coefs <- brain_growing@history$coefficients
  his.m <- brain_growing@history$mat.connectivity
  
  # colnames(his.coefs) <- c("Owner", "Seed",
  #                      "Round", "Epsilon Proportion",
  #                      "a Proportion", "Rewiring",
  #                      "Clustering", "Efficiency",
  #                      "Small World", "Modularity",
  #                      "Avg Path Length")
  
  time_start <- Sys.time()
  new.a <- now.a
  new.m <- now.m
  new.updates <- now.updates
  new.rewires <- now.rewires
  
  for(r in 1:n_rewires) {
    # updating nodes for n_update times
    for(u in 1:n_updates) {
      new.a <- trial_logistic(new.a, new.m,
                              eps = eps, a_param = a)
    }
    
    # incrementing age
    new.updates <- new.updates + n_updates
    new.rewires <- new.rewires + 1
    
    # rewiring
    new.m <- my_rewire(new.a,
                       new.m)
    # adding activities and clustering coefficient to the history
    his.a <- his.a %>% rbind(new.a)
    
    # saving the snapshot
    if(!(new.rewires %% freq_snapshot)){
      his.m[[new.rewires]] <- new.m
      new.coefs <- netmeas_coefs(brain_growing@initial$mat.connectivity[[1]],
                                 new.m,
                                 brain_growing@parameters,
                                 name = brain_growing@name,
                                 t_ = new.rewires)
      his.coefs <- his.coefs %>% rbind(new.coefs)
      }
    
    if(!quiet & !(new.rewires %% 100)) print(paste(brain_growing@name, "is now",
                           # new.updates, "updates and",
                           new.rewires, "rewires old."))
  }
  
  # updating the brain
  ## age
  brain_growing@age$updates <- new.updates
  brain_growing@age$rewires <- new.rewires
  
  ## history
  brain_growing@history$activities <- his.a
  brain_growing@history$mat.connectivity <- his.m
  brain_growing@history$coefficients <- his.coefs
  
  ## now
  brain_growing@now$activities <- new.a
  brain_growing@now$mat.connectivity <- new.m
  brain_growing@now$coefficients <- new.coefs
  
  (time_taken <- Sys.time() - time_start) %>% paste("for", brain_growing@name) %>%
    print()

  if(save_brain) save_vars("brain_growing", prefix = paste0("vat_", brain_growing@name))
  
  brain_growing %>% return()
  
}


# summary of aged brains --------------------------------------------------

trial_summary <- function(aged_brain){
  
  activities <- aged_brain@history$activities
  cl.c <- aged_brain@history$coef.clustering
  # cl.c %>% mean()
  # cl.c %>% plot()#main = paste("Clustering Coefficient of", b@name))
  
  cl.c_range <- cl.c %>% range()
  coef.clustering_normalized <- (cl.c_range[2]-cl.c_range[1])/mean(cl.c)
  
  # variance of activity of nodes at each update
  (variance_between_node <- apply(activities, 1, var)) # %>% plot(main = "between")
  
  # variance of each node at it's life time
  (variance_within_node <- apply(activities, 2, var)) # %>% plot(main = "within")
  
  # trial_summary <- list(name = aged_brain@name,
  #                       age = aged_brain@age,
  #                       parameters = aged_brain@parameters,
  #                       mean_variance = variance_between_node %>% mean(),
  #                       coef.clustering_normalized = coef.clustering_normalized
  #                       )
  
  output <- c(aged_brain@name,
              aged_brain@age$update,
              aged_brain@age$rewire,
              aged_brain@parameters$n_nodes,
              aged_brain@parameters$n_edges,
              aged_brain@parameters$eps,
              aged_brain@parameters$seed,
              variance_between_node %>% mean(),
              coef.clustering_normalized
  ) %>% t() %>% as.data.frame()
  
  colnames(output) <- c("name", "age_update", "age_rewire",
                        "n_nodes", "n_edges", "eps",
                        "seed", "mean_variance", "coef.clustering_normalized")
  
  return(output)
  
}


# make a code for brain ---------------------------------------------------

make_brain.code <- function(p_, name = NULL, b = NULL){
  # if(!is.null(b)) p_ <- b@parameters; name <- b@name
  
  r_ <- p_$round
  n_n <- p_$n_nodes/1000
  n_e <- p_$n_edges/1000
  
  # if(length(p_$params.eps_a)==6){
    e_ <- p_$params.eps_a[1:3] %>% paste(collapse = "v")
    a_ <- p_$params.eps_a[4:6] %>% paste(collapse = "v")
    
    seed <- paste0(# p_$params.eps_a[1:3],
                   # p_$params.eps_a[4:6],
                   2*n_n,
                   15*n_e,
                   100+r_,
                   collapse = "")
    seed <- gsub("[.]","",seed) %>% as.numeric()
    seed <- seed %% .Machine$integer.max
  # }
  # if(length(p_$params.eps_a)==4){
  #   e_ <- p_$params.eps_a[1:2] %>% paste(collapse = "v")
  #   a_ <- p_$params.eps_a[3:4] %>% paste(collapse = "v")
  #   
  #   seed <- paste0(p_$params.eps_a[1:2],
  #                  p_$params.eps_a[3:4],
  #                  #n_n,
  #                  #n_e,
  #                  r_,
  #                  collapse = "")
  #   seed <- gsub("[0.]","",seed) %>% as.numeric()
  #   seed <- seed %% .Machine$integer.max
  #   }
  
  
  
  # bc_ brain code
  bc_ <- paste0("eps-", e_,
                   "_a-", a_,
                   "_r-", r_)
  # number of nodes and edges
  n_e <- paste0("g-",
                n_n,
                # "k_",
                "k-",
                n_e,
                "k")
  # since the names should be different for eps_a parameters,
  # it needs a specific seed for itself.
  seed.name <- p_$params.eps_a %>%
    paste(collapse = "") %>%
    paste0(seed) %>% as.numeric()
  seed.name <- seed.name %% .Machine$integer.max
  list(name = give_name(seed = seed.name),
       braincode = paste(bc_, n_e, sep = "_"),
       seed = seed) %>% return()
  
}

# logistic update for trial -----------------------------------------------

trial_logistic <- function(a, m, eps, a_param = 1.7) {
  # unit.vector allows to calculate M_i by multiplying it the connectivity matrix
  if(length(eps)==1) eps <- rep(eps, length(a))
  eps <- eps %>% as.matrix()
  unit.vector <- matrix(1, length(a), 1)
  M <- m %*% unit.vector
  fx <- a %>% mini_logistic(a = a_param) %>% as.matrix() %>% t()
  (a_next <- (1 - eps)*fx + m %*% fx*eps / M) %>% t() %>% return()
  
}


