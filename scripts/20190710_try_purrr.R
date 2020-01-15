rm(list = ls())

setwd("I:/Thesis_Codes/Thesis_R_Codes/")

source("./functions/functions_extract.R")


path.to.brains <- "./data/5200-edges"
path.to.snp <- "./data/5200-edges-snapshots"


r_ <- 1

pat.tmp <- NULL

pat.tmp <- "_g-0.3k-5.2k"

pattern <- paste0("r-", r_, pat.tmp)
pattern <- pat.tmp

paste("Now reading round", r_, "of", pat.tmp) %>% print()


brain_files <- list.files(path = path.to.brains,
                          pattern = "*.RData")

r.this <- brain_files[grepl(pattern, brain_files)]
r.this <- r.this[grepl("_g-", r.this)]
r.this <- r.this[grepl(paste0("_r-", r_), r.this)]
r.this <- r.this[!grepl("coefs.wb.", r.this)]
brain_locations <- path.to.brains %>% 
  paste(r.this, sep = "/")

this.brain_location <- brain_locations#[1]

Sys.time() %>% print()
system.time(snp <- this.brain_location %>%
              ldply(extract_brains,
                    snapshots =  seq(5e3, 100e3, 5e3)
                    )
            ) %>% print()
save_vars(
  list.of.vars = "snp",
  prefix = paste("snp_round-", r_),
  path = path.to.snp
)






# load(this.brain_location)
# 
# l.m <- brain_case@history$mat.connectivity
# l.m[1] <- NULL
# l.m.compact <- l.m %>% compact()
# 
# b <- brain_case
# 
# 
# m <- l.m.compact[[500]]
# 
# 
# 
# 
# 
# 
# x <- list(m,m,m,m) %>% as_tribble()
# 
# xx <- x %>% nest() %>% mutate(m)
# 
# 
# 
# 




# library(purrr)
# 
# n_iris <- iris %>%
#   group_by(Species) %>%
#   nest()
# 
# 
# 
# mod_fun <- function(df) lm(Sepal.Length ~ ., data = df)
# 
# m_iris <- n_iris %>%
#   mutate(model = map(data, mod_fun))
# 

m <- brain_case@now$mat.connectivity
snapshot.interval <- 50e3
1e6/snapshot.interval
m.l <- list()
for(i in 1:(1e6*50*10/snapshot.interval)) m.l[[i]] <- m


