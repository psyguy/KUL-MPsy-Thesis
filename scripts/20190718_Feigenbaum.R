source("./functions/functions_extract.R")
library(tidyverse)
library(compiler) ## requires R >= 2.13.0
library(latex2exp)
library(tikzDevice)

logistic.map <- function(r, x, N, M){
  ## r: bifurcation parameter
  ## x: initial value
  ## N: Number of iteration
  ## M: Number of iteration points to be returned
  z <- 1:N
  z[1] <- x
  for(i in c(1:(N-1))){
    ## our own logistic map
    z[i+1] <- (1 - r *z[i]*z[i])
    ## the conventional logistic map
    # z[i+1] <- r *z[i]*(1 - z[i])
  }
  ## Return the last M iterations 
  z[c((N-M):N)]
}

logistic.map <- cmpfun(logistic.map) # same function as above

a.seq <- seq(0, 2.0,
             by=0.001)
# a.seq <- seq(1.6, 2.0,
#              by=0.0001)
N <- 4000
M <- 150
start.x <- 0.1

`Final Value` <- sapply(a.seq,
                        logistic.map,
                        x = start.x,
                        N = N,
                        M = M) %>% 
  as.vector()

a <- sort(rep(a.seq, (M+1)))
# plot(`Final Value` ~ r, pch=".", col=rgb(0,0,0,0.05))

d <- data.frame(a, `Final Value`)
Sys.time()
# system.time()

p <- ggplot(aes(a,
                `Final Value`),
            data = d) +
  # ylim(-1, 1) +
  # annotate("rect",
  #          xmin = 1.7,
  #          xmax = 1.9,
  #          ymin = -Inf,
  #          ymax = Inf,
  #          fill = "deepskyblue3",
  #          alpha = 0.3) +
  cowplot::theme_half_open() +
  geom_vline(xintercept = c(1.7, 1.8, 1.9),
             color = "orangered",
             size = 2,
             alpha = .9) +
  theme_bw() +
  theme(text = element_text(size = 25,
                            family = "CMU Classical Serif"),
        # axis.title.x = "\n $\\alpha$",
        # axis.title.y = "$x_t$",
        plot.title = element_blank(),
        legend.position = "none"
  ) +
  geom_point(size = 0.01,
             color = "azure4",
             alpha = 0.7) +
  scale_x_continuous(breaks = seq(0, 2, 0.5),
                   limits = c(0, 2)) +
  scale_y_continuous(breaks = seq(-1, 1, 0.5),
                     limits = c(-1, 1)) +
  coord_fixed() +
  ylab(TeX("$X_t$")) +
  xlab(TeX("$\n \\alpha$"))


cowplot::save_plot("Feigenbaum_0-to-2_2.png",
          p,
          base_height = 20,
          base_width = 20)

tikz(file = "plot_test.tex", width = 8.27, height = 11.69)

p
# print(plot.rc.all)

dev.off()


system.time(p)

ggsave("Feigenbaum_0-to-2_2.png",
       plot = p,
       width = 30,
       height = 30,
       dpi = 800,
       units = "cm")
Sys.time()

q <- 0.9^(1/15)

1-(0.9 + (15*(1-q)*q^14))

# plot(`Final Value`~a, type = ".")
