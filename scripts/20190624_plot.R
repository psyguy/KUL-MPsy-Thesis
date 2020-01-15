rm(list = ls())

source("./functions/functions_trial.R")
source("./functions/functions_netmeas.R")


# for loop over all brain cases -------------------------------------------


sampled.path <- "./data/"
sampled.names <- list.files(path = sampled.path, pattern = "*.RData")

coefs_all <- NULL
rewires <- seq(1, 10001, 5)
t <- Sys.time()
for(sampled in sampled.names){
  load(paste(sampled.path, sampled, sep = ""))
  brain_case@name %>% print()
  coefs_all <- coefs_all %>% rbind(brain_case %>% netmeas_coefs())
  
  plot.title<- paste(brain_case@name)

  sbs_h <- (brain_case@history$activities[rewires,1:300]) %>%
    as.data.frame()

  sbs_h[,12] %>% plot(main = brain_case@parameters$params.eps_a)

  data <- sbs_h %>% cbind(rewires)

  # Molten$variable %>% str()

  # data <- data.frame(time = seq(0, 23), noob = rnorm(24), plus = runif(24), extra = rpois(24, lambda = 1))
  Molten <- reshape::melt(data, id.vars = "rewires")

  ggplot(Molten,
         aes(x = rewires,
             y = value,
             colour = variable)) +
    geom_line(alpha = 0.5) + theme(legend.position="none") + ylim(-1, 1) + ggtitle(plot.title)
  # ggsave(paste0("activities.",brain_case@name,".png"), dpi = 500)
}
Sys.time()-t


# plotting coefficients over time -----------------------------------------


coefs_all %>% ggplot(aes(x = rewiring,
                         y = coef.clustering,
                         colour = a,
                         linetype = eps)) +
  geom_line(size = .75, alpha = 0.8) +  
              ggtitle("Clustering Coefficient")
ggsave("coef.clustering.png")

coefs_all %>% ggplot(aes(x = rewiring,
                         y = coef.efficiency,
                         colour = a,
                         linetype = eps)) +
  geom_line(size = .75, alpha = 0.8) + 
                       ggtitle("Global Efficiency")
ggsave("coef.efficiency.png")

coefs_all %>% ggplot(aes(x = rewiring,
                         y = coef.smallworld,
                         colour = a,
                         linetype = eps)) +
  geom_line(size = .75, alpha = 0.8) + 
  ggtitle("Small World index")
ggsave("coef.smallworld.png")

coefs_all %>% ggplot(aes(x = rewiring,
                         y = coef.modularity,
                         colour = a,
                         linetype = eps)) +
  geom_line(size = .75, alpha = 0.8) + 
  ggtitle("Modularity (fast greedy clustering)")
ggsave("coef.modularity.png")

coefs_all %>% ggplot(aes(x = rewiring,
                         y = coef.avgpathlength,
                         colour = a,
                         linetype = eps)) +
  geom_line(size = .75, alpha = 0.8) + 
  ggtitle("Average path length")
ggsave("coef.avgpathlength.png")


# saving the parameters ---------------------------------------------------

coefs_all$eps <- coefs_all$eps %>% as.character() %>% as.numeric()

parameters_df <- coefs_all %>%
  filter(rewiring == 200) %>%
  # arrange(eps) %>%
  unique() %>% 
  select(1:5)


# --- Graph 1 : If you want ONLY the table in your image :
# First I create an empty graph with absolutely nothing :
qplot(1:10, 1:10, geom = "blank") + theme_bw() + theme(line = element_blank(), text = element_blank()) +
  # Then I add my table :
  annotation_custom(grob = tableGrob(parameters_df))
ggsave("brain_parameters.png")

write.csv(parameters_df, file = "brain_parameters_hpc_20190625")

save_vars(list.of.vars = c("parameters_df", "coefs_all"), prefix = "9brains")


# plotting connectivity matrix --------------------------------------------

m_now <- brain_case@starting_values$mat.connectivity[[1]]
seriated <- m_now %>% seriate()
pimage(m_now, seriated)


# Implementing Ilias' code - didn't work :( -------------------------------

A <- m_now
deg <- m_now %>% rowSums()

# if there is a degree 0, in the inversion it stays 0
deg2 <- 1/sqrt(deg)
deg2[deg==0] <- 0


deginv <- deg2# %>% t() #%>% c(0)
L = diag(length(deg2)) - ((A %*% deginv) %*% t(deginv))  # Get the normalized Laplacian

# decompose the matrix to its eigenvectors/eigenvalues
# eigval, eigvec = np.linalg.eigh(L)

e <- L %>% eigen(symmetric = TRUE)

x <- e$values %>% sort(index.return=TRUE)
x$x %>% tail()

# takes second eigenvalue. The first is trivial solution. Next way to take more than one eigenvectors and do clustering
x <- e$vectors[2,] %>% sort(index.return=TRUE)
b <- A[x$ix,]
b <- b[,x$ix]
pimage(b)


# look at logistic --------------------------------------------------------

x <- seq(0,1, by=.01)
x^2 %>% plot()
x %>% mini_logistic(a=2.7) %>% plot(x=x,
                                       type = "l",
                                       xlim = c(0,1),
                                       ylim = c(-1,1))


# look at the distributions -----------------------------------------------

set.seed(1)

(x <- rbeta(200000, 1,1)) %>% densityplot()


bdis <- function(x, alpha, beta){
  
  x^(alpha-1)*(1-x)^(beta-1)/beta(alpha, beta)
  
}

alpha <- c(0.5, 1, 1, 1, 2, 5)# %>% rep(times = nrounds)
beta <- c(0.5, 1, 2, 5, 1, 1)# %>% rep(times = nrounds)
alphabeta <- data.frame(alpha,beta)

d <- matrix(0, ncol = 6, nrow = 1001) %>% data.frame()
for(i in 1:6){
  d[,i] <- seq(0,1,by=0.001) %>% bdis(alphabeta[i,1],alphabeta[i,2])
}
  

