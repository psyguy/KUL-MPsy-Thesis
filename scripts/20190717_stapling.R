#!/usr/bin/env Rscript

# Args <- commandArgs(TRUE)
# rounds <- as.numeric(Args[1])
rm(list = ls())

source("./functions/functions_extract.R")

pdf.names <- list.files(path = "./figures/",
                        pattern = "Final Rich")

pdf.names <- pdf.names[grepl("at 1000k rewirings_connectivities",
                             pdf.names)]
names <- gsub(" at 1000k rewirings_connectivities",
              "",
              pdf.names)
i <- 21
Sys.time()
for(i in seq_along(names)){
  print(i)
  if(!grepl("Jana Wouters", pdf.names[i])) next()
  pages <- (1+ (i-1)*40):(i*40)
  staplr::select_pages(pages,
                       "./figures/Network structure - all 25k.pdf",
                       paste("./figures/Network maturation of", names[i]))
}
Sys.time()



# making pdfs of only adj mats --------------------------------------------

pdf.names <- list.files(path = "./figures/connectivity-plots_25k",
                        pattern = "*.pdf")

pdf.names <- pdf.names[grepl("rewirings_connectivities",
                             pdf.names)]

names <- gsub(" at.*",
              "",
              pdf.names) %>% 
  unique() %>% 
  paste0(".pdf")

Sys.time()
for(i in seq_along(names)){
  print(i)
  pages <- (1+ (i-1)*40):(i*40)
  staplr::staple_pdf(input_files = paste0("./figures/connectivity-plots_25k/",
                                          pdf.names[pages]),
                     output_filepath = paste("./figures/Evolution of connectivities in",
                                             names[i])
                                             )
}
Sys.time()





path.to.pdfs <- "./figures/netstats-plots_25k"
coef.files <- list.files(path = path.to.pdfs, pattern = "Final Rich")

p <- coef.files[grepl("homo", coef.files)] %>% 
  c(coef.files[grepl("hypo-cha", coef.files)]) %>% 
  c(coef.files[grepl("hyper-cha", coef.files)]) %>% 
  c(coef.files[grepl("hypo-cou", coef.files)]) %>% 
  c(coef.files[grepl("hyper-cou", coef.files)])



staple_pdf(input_files = paste(path.to.pdfs, p, sep = "/"),
           output_filepath = NULL)

?savePlot
