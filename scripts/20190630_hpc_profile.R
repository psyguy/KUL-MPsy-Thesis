#!/usr/bin/env Rscript

Args <- commandArgs(TRUE)
index <- as.numeric(Args[1])

# rm(list=ls())

source("./functions/functions_reports.R")
options(bitmapType='cairo')
require(grid)

path.to.save <- "figures/20190630_hpc_figures/"
path.to.read <- "figures/20190629_hpc_figures/"
width.column.report <- 2240

# reading the aggregated harvested data
load("./data/coefs.all_20190630_1902.RData")

color.bg <- "white"
color.mino <- "deepskyblue3" #"blue2"
color.majo <- "orangered" #"brown3"
color.inter <- "olivedrab2" #"green4"
color.whole <- "azure4"

max.cl <- coefs.all$Clustering %>% max()
max.sw <- coefs.all$`Small World` %>% max()
max.mo <- coefs.all$Modularity %>% max()
max.de <- coefs.all$`Edge Density` %>% max()
max.pl <- coefs.all$`Avg Path Length` %>% max()
max.ef <- coefs.all$Efficiency %>% max()

names.list <- coefs.all$Owner %>% unique()

t3 <- Sys.time()
name <- names.list[index] %>% as.character()

print(name)
coefs.this.round <- coefs.all %>% filter(Owner == name)

title <- paste0(name,
                " (",
                coefs.this.round$`Verbal Description`[1],
                ")")

p.cl <- ggplot(data = coefs.this.round,
               aes(x = Rewiring,
                   y = Clustering,
                   colour = Partition)) +
  geom_line(size = .75, alpha = 0.7) +
  scale_colour_manual(values = c(color.inter, color.majo,
                                 color.mino, color.whole)) +
  ggplot2::ylim(0,max.cl) + 
  theme(legend.title = element_blank(),#element_text(size=12, color = "thistle4"),
        legend.justification=c(0,1), 
        legend.key.width = unit(1.5, "line"),
        legend.position= c(0, 1),#c(0.05, 0.95),
        legend.background = element_blank(),
        legend.text = element_text(size=8, color = "black", face="bold"),
        legend.key = element_blank())
# + labs(subtitle="Legend: Top-Left Inside the Plot")

p.sw <- ggplot(data = coefs.this.round,
               aes(x = Rewiring,
                   y = `Small World`,
                   colour = Partition)) +
  geom_line(size = .75, alpha = 0.7) +
  scale_colour_manual(values = c(color.inter, color.majo,
                                 color.mino, color.whole)) +
  ggplot2::ylim(0,max.sw) + 
  theme(legend.title = element_blank(),#element_text(size=12, color = "thistle4"),
        legend.justification=c(0,1), 
        legend.key.width = unit(1.5, "line"),
        legend.position= c(0, 1),#c(0.05, 0.95),
        legend.background = element_blank(),
        legend.text = element_text(size=8, color = "black", face="bold"),
        legend.key = element_blank())


p.mo <- ggplot(data = coefs.this.round,
               aes(x = Rewiring,
                   y = Modularity,
                   colour = Partition)) +
  geom_line(size = .75, alpha = 0.7) + 
  scale_colour_manual(values = c(color.inter, color.majo,
                                 color.mino, color.whole)) +
  ggplot2::ylim(0,max.mo) + 
  theme(legend.title = element_blank(),#element_text(size=12, color = "thistle4"),
        legend.justification=c(0,1), 
        legend.key.width = unit(1.5, "line"),
        legend.position= c(0, 1),#c(0.05, 0.95),
        legend.background = element_blank(),
        legend.text = element_text(size=8, color = "black", face="bold"),
        legend.key = element_blank())


p.ed <- ggplot(data = coefs.this.round,
               aes(x = Rewiring,
                   y = `Edge Density`,
                   colour = Partition)) +
  geom_line(size = .75, alpha = 0.7) + 
  scale_colour_manual(values = c(color.inter, color.majo,
                                 color.mino, color.whole)) +
  ggplot2::ylim(0,max.de) + 
  theme(legend.title = element_blank(),#element_text(size=12, color = "thistle4"),
        legend.justification=c(0,1), 
        legend.key.width = unit(1.5, "line"),
        legend.position= c(0, 1),#c(0.05, 0.95),
        legend.background = element_blank(),
        legend.text = element_text(size=8, color = "black", face="bold"),
        legend.key = element_blank())


p.pl <- ggplot(data = coefs.this.round,
               aes(x = Rewiring,
                   y = `Avg Path Length`,
                   colour = Partition)) +
  geom_line(size = .75, alpha = 0.7) + 
  scale_colour_manual(values = c(color.inter, color.majo,
                                 color.mino, color.whole)) +
  ggplot2::ylim(0,max.pl) + 
  theme(legend.title = element_blank(),#element_text(size=12, color = "thistle4"),
        legend.justification=c(0,1), 
        legend.key.width = unit(1.5, "line"),
        legend.position= c(0, 1),#c(0.05, 0.95),
        legend.background = element_blank(),
        legend.text = element_text(size=8, color = "black", face="bold"),
        legend.key = element_blank())


p.ef <- ggplot(data = coefs.this.round,
               aes(x = Rewiring,
                   y = Efficiency,
                   colour = Partition)) +
  geom_line(size = .75, alpha = 0.7) + 
  scale_colour_manual(values = c(color.inter, color.majo,
                                 color.mino, color.whole)) +
  ggplot2::ylim(0,max.ef) + 
  theme(legend.title = element_blank(),#element_text(size=12, color = "thistle4"),
        legend.justification=c(0,1), 
        legend.key.width = unit(1.5, "line"),
        legend.position= c(0, 1),#c(0.05, 0.95),
        legend.background = element_blank(),
        legend.text = element_text(size=8, color = "black", face="bold"),
        legend.key = element_blank())


figure <- ggarrange(p.cl, p.ed, p.sw, p.pl, p.mo, p.ef,
                    ncol=2, nrow=3,
                    common.legend = F
                    # legend="bottom"
)

require(grid)
top.text <- textGrob(paste("Evolution of coefficients for", title, "\n"),
                     gp=gpar(fontsize=25,font=8))
annotate_figure(figure, top = top.text)

paste0(path.to.save,"wb_", title, ".png") %>% ggsave(width = 14,
                                       height = 21,
                                       dpi = "retina")

paste("Plotting coefficients of",
      title,
      "took",
      (Sys.time()-t3)) %>% print()


# # reading network plots and gluing them to the coefficients over t --------
# 
# plot.names.nets <- list.files(path = path.to.read, pattern = "*rewirings.png")
# r.this <- plot.names.nets[grepl(name, plot.names.nets)] %>% sort()
# # to send the 1e+50 to the end
# r.this <- c(r.this[1], r.this[3:length(r.this)], r.this[2])
# 
# plot.name.coefs <- list.files(path = path.to.save, pattern = name)
# 
# plots.with.paths <- c(paste0(path.to.save, plot.name.coefs),
#                       paste0(path.to.read, r.this))
# 
# n.plots <- length(plots.with.paths)
# 
# t5 <- Sys.time()
# all.plots <- plots.with.paths %>%
#   lapply(png::readPNG) %>% 
#   lapply(grid::rasterGrob)
# 
# plot.n_col <- 3
# plot.n_row <- 2
# 
# Sys.time() - t5
# require(grid)
# paste0(path.to.save, "Profile of ", title, ".png") %>%
#   png(width = width.column.report*2*plot.n_col,
#       height = width.column.report*3*plot.n_row,
#       res = 400)
# gr.net <- gridExtra::grid.arrange(grobs=all.plots,
#                                   ncol = plot.n_col,
#                                   nrow = plot.n_row,
#                                   top = textGrob(paste(#"\n",
#                                                        "Profile of",
#                                                        title,
#                                                        "\n"),
#                                                  gp=gpar(fontsize=70,font=8)
#                                                  )
#                                   )
# dev.off()
# Sys.time() - t5
# 
# paste("Profile of",
#       title,
#       "took",
#       (Sys.time()-t5)) %>% print()
# 
# 