#!/usr/bin/env Rscript

# Args <- commandArgs(TRUE)
# rounds <- as.numeric(Args[1])
rm(list = ls())

source("./functions/functions_extract.R")

path.to.pdfs <- "./figures/pdfs_1-2/"

con.files <- list.files(path = path.to.pdfs, pattern = "*.pdf")

i <- 7

this.con.files <- con.files[(3*(i-1)+1):(3*i)]



load("./data/snp.lean_all_5k_20190713_1356.RData")
snp <- snp.lean
rm(snp.lean)

all.owners <- snp$Owner %>% 
  unique() %>% 
  sort() %>% 
  as.character() %>%
  as.list()

l <- paste0(path.to.pdfs, this.con.files) #%>% 
as.list() %>% 
  map(~image_read_pdf(.x))


whole2 <- image_append(c(image_read_pdf(l[1]),
                         image_read_pdf(l[3])))

plot(image_read_pdf(l[1]))
plot(whole2)
savePlot("whole2.pdf", type = "pdf") 

panel.upr <- image_read_pdf("whole2.pdf")
final <-  image_append(c(image_append(panel.upr),l[[1]]), stack = TRUE)
plot(final)
savePlot("final.pdf", type = "pdf") 





whole <- c(l[[3]], l[[2]]) %>% 
  image_append() %>% 
  c(l[[1]]) %>% 
  image_append(stack = TRUE)

plot(whole)
savePlot(type = "pdf")

graph2pdf(width = 10, height = 5)

library(pdftools)
pdftools::pdf_combine()

library(grImport)
grImport

staple_pdf("./figures/Network statistics by 25k - 5200 edges")
