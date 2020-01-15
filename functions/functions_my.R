# Intro -------------------------------------------------------------------

## This script keeps all functions needed for GongvLeeuwen2004.
## It's standalone, doesn't need to call any other script

# loading packages --------------------------------------------------------

list.of.packages <- c("tidyverse",
                      "dplyr",
                      "plyr",
#                      "corrplot",
                      "Hmisc",
		      "png",
                      "seriation",
                      "ggpubr",
		      "HHG",
                      "magick",
                      "export",
                      "staplr",
                      "igraph",
                      "ggplot2",
                      "gplots",
                      "reshape",
                      "brainGraph",
                      "gridExtra")
new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) {
  install.packages(new.packages)#, repos='https://lib.ugent.be/CRAN/')
}
tmp <- lapply(list.of.packages, require, character.only = TRUE)

rm(list.of.packages, new.packages, tmp)


# mini functions ----------------------------------------------------------

mini_logistic <- function(x, a = 1.7) {
  1 - a * (x ^ 2)
}

my_gsub <- function(str, pat, sub = "") gsub(pat, sub, str)



# vectorizing symmetrical adj. matrix and vice versa ----------------------


vec2mat <- function(v_){
  if(is.list(v_)) v_ <- v_[[1]]
  l_ <- length(v_)
  n_ <- (1 + sqrt(1 + 8*l_))/2
  m_ <- matrix(0, n_, n_)
  m_[lower.tri(m_, diag = FALSE)] <- v_
  (m_ + t(m_)) %>% return()
}


mat2vec <- function(m_){
  l_ <- nrow(m_)
  m_[m_ %>% lower.tri(diag = FALSE)] %>%
    as.vector() %>%
    return()
}
# saves a backup of variables ---------------------------------------------

save_vars <- function(list.of.vars = NULL,
                      envir = parent.frame(),
                      prefix = "StatusQuo",
                      path = "data") {
  mine_gsub <- function(x, pat, rep)
    gsub(pat, rep, x)
  
  if(is.null(list.of.vars)) list.of.vars <- ls(envir = envir) # setdiff(ls(), lsf.str())
  
  date_time <- Sys.time() %>%
    mine_gsub(" ", "_") %>%
    mine_gsub("-", "") %>%
    mine_gsub(":", "") %>%
    strtrim(13) # getting rid of timezone, etc.
  
  if (!is.null(path))
    path <- paste0(path, "/")
  file_name <- paste0(path, prefix, "_", date_time, ".RData")
  
  save(list = list.of.vars, file = file_name)
}


# to make random graphs of certain size and #links -------------------------

make_random_graph <- function(size = 5,
                              num.links = 8,
                              distribution = "binary",
                              parameters = c(0, 1),
                              seed = rnorm(1)) {
  set.seed(seed)
  edges <- rep(0, size * (size - 1) / 2)
  
  if (distribution == "binary") {
    ones <- sample.int(length(edges), num.links)
    edges[ones] <- 1
  }
  
  g <- matrix(0, size, size)
  g[base::lower.tri(g, diag = FALSE)] <- edges
  
  (adj.matrix <- g + t(g)) %>% return()
}


# swapping the desired edges ----------------------------------------------

swap_edges <- function(connectivity.matrix,
                       from.1,
                       to.1,
                       from.2,
                       to.2) {
  edge.1 <- connectivity.matrix[from.1, to.1]
  edge.2 <- connectivity.matrix[from.2, to.2]
  
  connectivity.matrix[from.1, to.1] <- edge.2
  connectivity.matrix[to.1, from.1] <- edge.2
  
  connectivity.matrix[from.2, to.2] <- edge.1
  connectivity.matrix[to.2, from.2] <- edge.1
  
  connectivity.matrix %>% return()
  
}

# my small rewiring function ----------------------------------------------

my_rewire <- function(a, m, global_minmax = FALSE, blind_swap = FALSE) {
  distances <- a %>% tail(1) %>% my_coherenceD()
  
  i_1 <- sample.int(ncol(m), 1) -> i_2
  d_ <- distances[, i_1]
  m_ <- m[, i_1]
  
  m_d_ <- (1 - m_)*d_
  m_d_[m_d_==0|is.na(m_d_)] <- 99
  
  j_connected_maxdist <- (m_*d_) %>% which.max() -> j_1
  j_disconnected_mindist <- m_d_ %>% which.min() -> j_2
  
  ## my wrong Way of doing it before 2019-06-19
  # j_1 <- which.min(d_)
  # j_2 <- which.max(d_)
  
  
  if(global_minmax){
    ind_min <- which(distances == min(distances, na.rm=T),
                     arr.ind = TRUE)[1,]
    ind_max <- which(distances == max(distances, na.rm=T),
                     arr.ind = TRUE)[1,]
    
    i_1 <- ind_min[1]
    j_1 <- ind_min[2]
    
    i_2 <- ind_max[1]
    j_2 <- ind_max[2]
  }
  
  if(blind_swap){
    m <- m %>% swap_edges(
      from.1 = i_1,
      to.1 = j_1,
      from.2 = i_2,
      to.2 = j_2)
    m %>% return()
  }
  
  m[i_1,j_1] <- m[j_1,i_1] <- 0
  m[i_2,j_2] <- m[j_2,i_2] <- 1
  
  m %>% return()
  
}


# coherence D -------------------------------------------------------------

my_coherenceD <- function(x.input) {
  x.duplicated <- rep(1, ncol(x.input)) %*% x.input
  d <- (t(x.duplicated) - x.duplicated) %>% abs()
  diag(d) <- NA
  d %>% return()
  
}


# make parameter distributions --------------------------------------------

make_paramdist <- function(alpha_beta,
                           range_param,
                           l = 300,
                           seed = -99){
  alpha_beta <- alpha_beta %>% as.matrix()
  ## if the alpha_beta parameter has a length of 3, it makes descreet distribution
  # if(length(alpha_beta)>2 | ncol(alpha_beta)>2){
  s_ <- alpha_beta %>% sum()
  n_low <- (as.numeric(alpha_beta[1])*l/s_) %>% round(digits=0)
  n_medium <- (as.numeric(alpha_beta[2])*l/s_) %>% round(digits=0)
  n_high <- (as.numeric(alpha_beta[3])*l/s_) %>% round(digits=0)
  
  v_low <- range_param[1] %>% rep(n_low)
  v_medium <- range_param %>% mean() %>% rep(n_medium)
  v_high <- range_param[2] %>% rep(n_high)
  
  v_ <- c(v_low, v_high, v_medium)
  v_ %>% return()
  # }
  # if(length(alpha_beta)<=2 | ncol(alpha_beta)<=2){
  #   alpha_beta <- alpha_beta %>% as.numeric()
  #   range_param <- range_param %>% sort()
  #   set.seed(seed)
  #   v_ <- l %>% rbeta(alpha_beta[1], alpha_beta[2])
  #   v_ <- v_*(range_param[2]-range_param[1]) + (range_param[1])
  #   v_ %>% return()
  # }
}


# make parameter vector ---------------------------------------------------

make_paramvect <- function(s = "0.3x1, 0.4x0, 0.5x0",
                           n = 300,
                           seed = -99){
  s <- gsub(" ", "", s)
  p_ <-  (s %>% strsplit(","))[[1]]
  
  v_ <- c()
  for(p in p_){
    l_r_ <- strsplit(p, "x")[[1]] %>% as.numeric()
    v_ <- v_ %>% c(rep(l_r_[1], l_r_[2]*n))
  }
  set.seed(seed)
  sample(v_) %>% return()
}


# clustering coeeficient (transitivity) -----------------------------------

my_clustceof <- function(m) {
  tr <- function(m)
    sum(diag(m), na.rm = TRUE)
  
  m <- m %>% as.matrix()
  if (dim(m)[1] != dim(m)[2])
    stop("m must be a square matrix")
  
  m2 <- m %*% m
  m3 <- m2 %*% m
  
  tr(m3) / (sum(m2) - tr(m2))
}


give_name <- function(num = 1, seed = -99, gender = "b"){
  
  set.seed(seed)
  
  firstnames_male <- c("Adrien","Alexander","Alexandre","Alexis","Anthony","Antoine",
                       "Arnaud","Arne","Arno","Arthur","Axel","Benjamin","Bram","Brent",
                       "Bryan","Clement","Corentin","Cyril","Cedric","Daan","David",
                       "Dorian","Dries","Dylan","Elias","Florian","Gilles","Guillaume",
                       "Hugo","Jarne","Jason","Jasper","Jens","Jonas","Jonathan","Jordan",
                       "Jordy","Julien","Justin","Jeremy","Kevin","Kobe","Lars","Lennert",
                       "Liam","Logan","Louis","Loic","Luca","Lucas","Lukas","Maarten",
                       "Martin","Mathias","Mathieu","Matthias","Maxim","Maxime","Mehdi",
                       "Michiel","Milan","Mohamed","Nathan","Nick","Nicolas","Niels","Noah",
                       "Olivier","Pierre","Pieter","Quentin","Quinten","Robbe","Robin",
                       "Romain","Ruben","Ryan","Sam","Samuel","Sander","Seppe","Simon",
                       "Stef","Stijn","Sebastien","Thibault","Thibaut","Thomas","Theo",
                       "Tibo","Tim","Tom","Tristan","Valentin","Victor","Vincent","Ward",
                       "William","Wout","Yannick")
  
  firstnames_female <- c("Alexandra","Alexia","Alice","Alicia","Aline","Amandine",
                         "Amber","Amelie","Anais","Anke","Anna","Anouk","Audrey",
                         "Aurelie","Axelle","Bo","Britt","Camille","Caro","Caroline",
                         "Charlotte","Chiara","Chloe","Clara","Celia","Celine",
                         "Delphine","Eline","Elisa","Elise","Ellen","Elodie","Emilie",
                         "Emma","Estelle","Eva","Fanny","Febe","Femke","Fien","Fiona",
                         "Fleur","Florence","Hannah","Hanne","Imane","Ine","Ines",
                         "Ines","Jade","Jana","Jessica","Jolien","Julie","Juliette",
                         "Justine","Kaat","Kato","Kelly","Lara","Laura","Laure",
                         "Lauren","Lien","Lies","Lisa","Lise","Lore","Lotte","Louise",
                         "Lucie","Luna","Lea","Manon","Margaux","Margot","Marie",
                         "Marine","Marthe","Melissa","Morgane","Melanie","Nina","Noa",
                         "Noemie","Oceane","Ophelie","Pauline","Rania","Sara","Sarah",
                         "Silke","Sofie","Sophie","Valentine","Victoria","Yana",
                         "Yasmine","Zoe")
  
  lastnames <- c("Adam","Aerts","Baert","Bauwens","Beckers","Bertrand","Bogaert",
                 "Bogaerts","Bosmans","Carlier","Charlier","Christiaens","Claes",
                 "Claessens","Claeys","Cools","Coppens","Cornelis","De Backer",
                 "De Clercq","De Cock","De Meyer","De Pauw","De Ridder","De Smedt",
                 "De Smet","De Vos","De Wilde","Declercq","Denis","Deprez","Desmet",
                 "Devos","Dubois","Dumont","Dupont","Evrard","Fontaine","Francois",
                 "Geerts","Goossens","Gerard","Hendrickx","Hermans","Jacobs","Jansen",
                 "Janssen","Janssens","Lambert","Lambrechts","Laurent","Lauwers",
                 "Leclercq","Lejeune","Lemaire","Lemmens","Lenaerts","Leroy","Maes",
                 "Martens","Martin","Mathieu","Mertens","Michel","Michiels","Moens",
                 "Noel","Pauwels","Peeters","Petit","Pieters","Renard","Segers",
                 "Simon","Simons","Smet","Smets","Stevens","Thomas","Thys",
                 "Timmermans","Van Damme","Van De Velde","Van Den Broeck","Van Dyck",
                 "Vandenberghe","Verbeke","Verheyen","Verhoeven","Verlinden",
                 "Vermeersch","Vermeiren","Vermeulen","Verschueren","Verstraete",
                 "Verstraeten","Wauters","Willems","Wouters","Wuyts")
  
  f_n <- firstnames_male %>% c(firstnames_female)
  if(gender=="m") f_n <- firstnames_male
  if(gender=="f") f_n <- firstnames_female
  
  f_n <- f_n %>% sample(num)
  l_n <- lastnames %>% sample(num)
  
  paste(f_n, l_n) %>% return()
  
}


