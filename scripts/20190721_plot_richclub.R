#!/usr/bin/env Rscript

# Args <- commandArgs(TRUE)
# rounds <- as.numeric(Args[1])
# rm(list = ls())

source("./functions/functions_extract.R")
load("./data/snp_only-1e+6_20190723_1404.RData")

# load("./data/snp.lean_all_5k_20190713_1356.RData")
# snp <- snp.lean
# rm(snp.lean)

all.owners <- snp$Owner %>% 
  unique() %>% 
  sort() %>% 
  as.character() %>%
  as.list()

all.partitions <- snp$Partition %>% 
  unique() %>% 
  sort() %>% 
  as.character() %>%
  as.list()

all.vd <- snp$Verbal.Description %>% 
  unique() %>% 
  sort() %>% 
  as.character() %>%
  as.list()

name.this.owner <- all.owners[[1]]

# rand.graphs <- snp$adj.mat.vect[[1]] %>% 
#   vec2mat() %>% 
#   graph_from_adjacency_matrix() %>% 
#   sim.rand.graph.par(200)

# The problematic owners (those who died out)
# are cause error in rendom graph generation,
# so the normalized Rich Club cannot be calculated for them.
# The indices are:
# 14 ("David De Ridder"), 17 ("Daan Dubois"),
# 41 ("Audrey Claeys"), 44 ("Chiara Bosmans")

t <- Sys.time()
system.time(all.owners[45:50] %>%
              map(extract_plot_rc.btwn,
                  snp = snp)
)
Sys.time() - t
# system.time(all.partitions %>%
#               map2(all.vd, extract_plotcoefs.glued,
#                   snp = snp)
# )

system.time(
  for(vd in all.vd){
    for(part in all.partitions){
      extract_plotcoefs.glued(this.Verbal.Description = vd,
                              this.Partition = part, 
                              snp=snp)
    }
  }
)




path.to.pdfs <- "./figures"
coef.files <- list.files(path = path.to.pdfs, pattern = "*.pdf")
coef.files <- coef.files[grepl("Rich Club", coef.files)]

p <- coef.files[grepl("homo", coef.files)] %>% 
  c(coef.files[grepl("hypo-cha", coef.files)]) %>% 
  c(coef.files[grepl("hyper-cha", coef.files)]) %>% 
  c(coef.files[grepl("hypo-cou", coef.files)]) %>% 
  c(coef.files[grepl("hyper-cou", coef.files)])



staple_pdf(input_files = paste(path.to.pdfs, p, sep = "/"),
           output_filepath = "partitioned per family.pdf")


?savePlot


# extract_plotcoefs.glued(this.Verbal.Description = all.vd[[3]],
#                         this.Partition = all.partitions[[3]], 
#                         snp=snp)
