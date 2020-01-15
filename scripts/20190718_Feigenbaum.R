source("./functions/functions_extract.R")

library(compiler) ## requires R >= 2.13.0
library(latex2exp)

logistic.map <- function(r, x, N, M){
  ## r: bifurcation parameter
  ## x: initial value
  ## N: Number of iteration
  ## M: Number of iteration points to be returned
  z <- 1:N
  z[1] <- x
  for(i in c(1:(N-1))){
    z[i+1] <- (1 - r *z[i]*z[i])
  }
  ## Return the last M iterations 
  z[c((N-M):N)]
}

logistic.map <- cmpfun(logistic.map) # same function as above

a.seq <- seq(0, 2.0,
             by=0.0005)
# a.seq <- seq(1.6, 2.0,
#              by=0.0001)
N <- 4000
M <- 200
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
  geom_point(shape = ".")

pp <- p + 
  ylim(-1, 1) +
  annotate("rect",
           xmin = 1.7,
           xmax = 1.9,
           ymin = -Inf,
           ymax = Inf,
           fill = "lightsteelblue1",
           alpha = .3) +
  
  geom_vline(xintercept = c(1.7, 1.8, 1.9),
             color = "orangered",
             alpha = .9,
             linetype = "longdash") +
  
  ylab(TeX("$x_t$")) +
  xlab(TeX("\n $\\alpha$"))

system.time(p)

ggsave("Feigenbaum_0-to-2.png",
       plot = pp,
       width = 30,
       height = 20,
       dpi = 800,
       units = "cm")
Sys.time()


# plot(`Final Value`~a, type = ".")
