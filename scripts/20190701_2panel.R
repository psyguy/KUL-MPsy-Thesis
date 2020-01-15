# rm(list=ls())
library(ggpubr)
library(grid)

source("./functions/functions_reports.R")

# path.to.coefs <- "figures/20190630_hpc_figures/"
path.to.pngs <- "figures/5400/"
width.panel <- 4480
height.panel <- 6720


names.coefs <- list.files(path = path.to.pngs, pattern = "wb_") %>% sort()
names.netvizs <- list.files(path = path.to.pngs, pattern = "after") %>% sort()


titles <- gsub("wb_","", names.coefs)
titles <- gsub(".png","", titles)

# i <- 1

# reading network plots and gluing them to the coefficients over t --------

t6 <- Sys.time()
for(i in 1:50){
  
  print(i)
  
  t5 <- Sys.time()
  plots.with.paths <- paste0(path.to.pngs,
                            c(names.coefs[i], names.netvizs[i]))
  two.pannels <- plots.with.paths %>%
    lapply(png::readPNG) %>%
    lapply(grid::rasterGrob)
  Sys.time() - t5
  
  require(grid)
  paste0(path.to.pngs, "2-panel profile of ", titles[i], ".png") %>%
    png(width = width.panel*2,
        height = height.panel,
        res = 400)
  gr.net <- gridExtra::grid.arrange(grobs=two.pannels,
                                    ncol = 2,
                                    top = textGrob(paste(#"\n",
                                                         "Profile of",
                                                         titles[i],
                                                         "\n"),
                                                   gp=gpar(fontsize=40,font=8)
                                                   )
                                    )
  dev.off()

  paste("2-panel profile of",
        titles[i],
        "took",
        (Sys.time()-t5)) %>% print()

}
Sys.time() - t6
