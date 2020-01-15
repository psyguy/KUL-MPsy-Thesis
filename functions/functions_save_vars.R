library(tidyverse)
save_vars <- function(list.of.vars = NULL,
                      envir = parent.frame(),
                      prefix = "StatusQuo",
                      path = "data") {
  mine_gsub <- function(x, pat, rep)
    gsub(pat, rep, x)
  
  if(is.null(list.of.vars)) list.of.vars <- ls(envir = envir) # setdiff(ls(), lsf.str())
  
  date_time <- Sys.time() %>%
    mine_gsub(" ", "_") %>%
    mine_gsub("-", "") %>%
    mine_gsub(":", "") %>%
    strtrim(13) # getting rid of timezone, etc.
  
  if (!is.null(path))
    path <- paste0(path, "/")
  file_name <- paste0(path, prefix, "_", date_time, ".RData")
  
  save(list = list.of.vars, file = file_name)
}