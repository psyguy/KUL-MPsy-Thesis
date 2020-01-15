# Intro -------------------------------------------------------------------

## This script keeps all functions needed for GongvLeeuwen2004.
## It's standalone and loads functions_my automatically

# rm(list = ls())


# loading the basic functions and packages --------------------------------

if(!exists("path.to.functions_my")) path.to.functions_my <- "functions"
source(paste0(path.to.functions_my,"/functions_my.R"))


# logistic function -------------------------------------------------------

GongvLeeuwen2004_logistic <- function(x.input,
                                      connectivity.matrix,
                                      a = 1.7,
                                      eps = 0.8) {
  # eps: coupling strength
  # eps <- 0.8
  # x.input <- x.init
  # connectivity.matrix <- conn
  
  x.now <- x.input <- x.init
  n.nodes <- x.now %>% ncol()
  
  
  M <- connectivity.matrix %*% rep(1, n.nodes) %>% t()
  fx <- x.now %>% mini_logistic() %>% as.matrix()# %>% t()
  
  neighbors <- fx %*% connectivity.matrix * eps / M
  neighbors[is.nan(neighbors)] <- 0
  
  x.next <- (1 - eps) * fx + neighbors
  
  x.next %>% return()
  
}


# calculating coherence D -------------------------------------------------

GongvLeeuwen2004_coherenceD <- function(x.input) {
  x.duplicated <- rep(1, ncol(x.input)) %*% x.input
  d <- (t(x.duplicated) - x.duplicated) %>% abs()
  diag(d) <- NA
  d %>% return()
  
}


# GongvLeeuwen2004 adaptive rewiring  ------------------------------------------------------

GongvLeeuwen2004_adaptive_rewire <- function(input) {
  x.start <- input$x.tot %>% tail(1)
  conn.start <- input$connectivity.matrix
  
  x.end <- x.start %>% GongvLeeuwen2004_logistic(conn.start)
  conn.end <- my_rewire(x.end, conn.start)
  
  g <-
    igraph::graph_from_adjacency_matrix(conn.end, mode = "undirected")
  ClCoef <- igraph::transitivity(g)
  
  output <- list(
    x.tot = rbind(input$x.tot, x.end),
    clustering.coefficient = c(input$clustering.coefficient, ClCoef),
    connectivity.matrix = conn.end
  )
  class(output) <- append(class(output), "GongvLeeuwen2004")
  output %>% return()
}
