# Intro -------------------------------------------------------------------

## This script keeps all functions needed for GongvLeeuwen2004.
## It's standalone and loads functions_my automatically

# rm(list = ls())


# loading the basic functions and packages --------------------------------

if(!exists("path.to.functions_my")) path.to.functions_my <- "functions"
source(paste0(path.to.functions_my,"/functions_my.R"))


# logistic function -------------------------------------------------------
# 
# heartbeat_logistic <- function(activities_,
#                                m, #connectiviti matrix
#                                a = 1.7,
#                                eps = 0.8) {
#   # eps: coupling strength
#   # eps <- 0.8
#   # x.input <- x.init
#   # connectivity.matrix <- conn
#   
#   activities_now <- activities_
#   num_nodes <- activities_now %>% ncol()
#   
#   
#   M <- m %*% rep(1, num_nodes) %>% t()
#   fx <- activities_now %>% mini_logistic() %>% as.matrix()
#   
#   neighbors <- fx %*% m * eps / M
#   neighbors[is.nan(neighbors)] <- 0
#   
#   activities_next <- (1 - eps) * fx + neighbors
#   
#   activities_next %>% return()
#   
# }


# making a S4 class -------------------------------------------------------

brain <- setClass("brain",
                    slots = list(name = "character",
                                 birthday = "character",
                                 age = "list",
                                 parameters = "list", # num_nodes, num_edges
                                 initial = "list", # activities
                                 history = "list", # activities, coef.clustering, connmatrix
                                 now = "list" # activities, coef.clustering, mat.connectivity
                                 ))

# setMethod("brain", signature(object = "brain"), function(object) {
#   object@sides
# })


# updating brain activities -----------------------------------------------

update <- function(object) {
  if(!testClass(object,"brain")) stop("Gimme a brain!")
  
  beat <- object@age$beat
  m <- object@now$mat.connectivity
  a <- object@now$activities
  c <- object@now$coef.clustering

  object@age$beat <- beat + 1  
  
  a_next <- object %>% heartbeat_logistic()
  
  object@history$activities <-
    object@history$activities %>%
    rbind(a_next)
  
  object@now$activities <- a_next

  
  object %>% return()
  }


# rewiring the brain ------------------------------------------------------

rewire <- function(object) {
  if(!testClass(object,"brain")) stop("Gimme a brain!")
  
  minute <- object@age$minute
  l.m <- object@now$mat.connectivity -> m
  if(is.list(l.m)) m <- l.m[[length(l.m)]]
  a <- object@now$activities
  c <- object@now$coef.clustering
  
  object@age$minute <- minute + 1
  
  m_next <- heartbeat_rewire(a, m)
  c_next <- m %>% my_clustceof()
  
  object@history$mat.connectivity[object@age$minute] <- m_next
  object@history$coef.clustering <- c(object@history$coef.clustering, c_next)
  
  object@now$mat.connectivity <- m_next
  object@now$coef.clustering <- c_next

  
  object %>% return()
  }


# logistic function -------------------------------------------------------

heartbeat_logistic <- function(b,
                               #a = 1.7,
                               eps = 0.8) {
  # eps: coupling strength
  n.node <- b@parameters$num_nodes
  l.m <- b@now$mat.connectivity -> m

  if(is.list(l.m)) m <- l.m[[length(l.m)]]
  a <- b@now$activities
  # unit.vector allows to calculate M_i by multiplying it the connectivity matrix
  unit.vector <- matrix(1, n.node, 1)
  M <- m %*% unit.vector
  fx <- a %>% mini_logistic() %>% as.matrix() %>% t()
  
  a_next <- (1 - eps)*fx + m %*% fx*eps / M# %>% t()
  
  a_next %>% t() %>% return()
  
}


# GongvLeeuwen2004 adaptive rewiring  ------------------------------------------------------


heartbeat_rewire <- function(a, m) {

  distances <- a %>% my_coherenceD()
  
  i_ <- sample.int(ncol(m), 1)
  d_ <- distances[, i_]
  
  j_1 <- which.min(d_)
  j_2 <- which.max(d_)
  
  # if(!conn[i_,j_1]) conn %>% return()
  m <- m %>% swap_edges(
    from.1 = i_,
    to.1 = j_1,
    from.2 = i_,
    to.2 = j_2
  )

  m %>% return()
}
  