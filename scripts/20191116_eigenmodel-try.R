install.packages("eigenmodel")
library("eigenmodel")

data(YX_Friend)

fit<-eigenmodel_mcmc(Y=YX_Friend$Y,X=YX_Friend$X,R=2,S=50,burn=50)

# in general you should run the Markov chain longer than 50 scans

plot(fit)


library(shiny)
runUrl("http://www.stat.washington.edu/hoff/sParXiv/sParXiv.zip")

install.packages("GraphAlignment")

install.packages("gtools")
