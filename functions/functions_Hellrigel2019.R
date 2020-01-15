# Intro -------------------------------------------------------------------

## This script keeps all functions needed for Hellrigel2019.
## It's standalone and loads functions_my automatically

# rm(list = ls())


# loading the basic functions and packages --------------------------------

if (!exists("path.to.functions_my"))
  path.to.functions_my <- "functions"
source(paste0(path.to.functions_my, "/functions_my.R"))


# logistic function -------------------------------------------------------


Hellrigel2019_logistic <- function(inp.x,
                                   connectivity.matrix,
                                   a = 1.7,
                                   eps = 0.8,
                                   order = 0) {
  inp.dim <- inp.x %>% dim()
  height <- inp.dim[1]
  leng <- inp.dim[2]
  if (!is.null(leng))
    x.now <- inp.x[, leng]
  if (is.null(height))
    height <- inp.x %>% length()
  if (is.null(leng))
    x.now <- inp.x
  
  # Higher orders are not implemented yet
  #if(order) x <- inp.x[,len-order:len]
  
  # unit.vector allows to calculate M_i by multiplying it the connectivity matrix
  unit.vector <- matrix(1, height, 1)
  M <- connectivity.matrix %*% unit.vector
  fx <- x.now %>% mini_logistic(a = a) %>% as.matrix()
  
  x.next <- (1 - eps) * fx + connectivity.matrix %*% fx * eps / M
  
  x.next %>% return()
  
}


# calculating csync error -------------------------------------------------

Hellrigel2019.syncerror <- function(x.input, connectivity.matrix) {
  inp.dim <- x.input %>% dim()
  height <- inp.dim[1]
  leng <- inp.dim[2]
  if (!is.null(leng))
    x.now <- x.input[, leng]
  if (is.null(height))
    height <- x.input %>% length()
  if (is.null(leng))
    x.now <- x.input
  
  x.duplicated <- x.now
  for (i in 1:(height - 1)) {
    x.duplicated <- x.now %>% cbind(x.duplicated)
  }
  
  x.duplicated <- x.duplicated %>% as.matrix()
  
  d <- (t(x.duplicated) - x.duplicated) %>% abs()
  
  diag(d) <- 0
  
  d %>% return()
  
}
