load("./data-pc/5200-edges/life-10_eps-0v6v0_a-1v5v0_r-7_g-0.3k-5.2k_Ruben Declercq_20190709_2253.RData")
load("./data-pc/5200-edges/life-10_eps-1v5v0_a-0v6v0_r-8_g-0.3k-5.2k_Margaux Bogaert_20190709_2221.RData")
load("./data-pc/5200-edges/life-10_eps-1v5v0_a-0v6v0_r-9_g-0.3k-5.2k_Hannah Bauwens_20190709_2225.RData")


# rm(list = ls())
load("./data-pc/5200-edges/life-10_eps-1v5v0_a-0v6v0_r-9_g-0.3k-5.2k_Hannah Bauwens_20190709_2225.RData")

# setwd("I:/Thesis_Codes/Thesis_R_Codes/")

source("./functions/functions_extract.R")


path.to.brains <- "./data-pc/5200-edges"
path.to.snp <- "./data/5200-edges-snapshots"
pattern <- "_g-0.3k-5.2k"
brain_files <- list.files(path = path.to.brains,
                          pattern = "*.RData")
brain_files_1e6 <- brain_files[grepl("life-10", brain_files)]

r.this <- brain_files[grepl("life-10", brain_files)]
r.this <- r.this[!grepl("coefs.wb.", r.this)]
brain_locations <- path.to.brains %>%
  paste(r.this, sep = "/")
path.to.snp <- "./data-pc/1e6-snapshots"
brain_files <- list.files(path = path.to.brains, pattern = "*.RData")
r.this <- brain_files[grepl("life-10", brain_files)]
brain_locations <- path.to.brains %>%
  paste(r.this, sep = "/")
b <- brain_case


for(fn in brain_files_1e6){
  
}


extract_now <- function(b){
  iden <- extract_identifiers(b)$whole
  owner <- iden$Owner
  verbal <- iden$Verbal.Description
  
  m <- b@now$mat.connectivity
  a <- b@now$activities
  
  o <- list(
    owner = owner,
    verbal = verbal,
    m = m,
    a = a
  )
  o %>% return()
}


dist_maker <- function(l.now){
  num.persons <- 50
  dist.m <- matrix(0, num.persons, nrow = num.persons)
  rownames(dist.m) <- colnames(dist.m) <- c(
    "Hyper-coupled minority", "Hyper-coupled minority",
    "Hyper-coupled minority", "Hyper-coupled minority",
    "Hyper-coupled minority", "Hyper-coupled minority",
    "Hyper-coupled minority", "Hyper-coupled minority",
    "Hyper-coupled minority", "Hyper-coupled minority",
    "Hyper-chaotic minority", "Hyper-chaotic minority",
    "Hyper-chaotic minority", "Hyper-chaotic minority",
    "Hyper-chaotic minority", "Hyper-chaotic minority",
    "Hyper-chaotic minority", "Hyper-chaotic minority",
    "Hyper-chaotic minority", "Hyper-chaotic minority",
    "Homogeneous Society", "Homogeneous Society",
    "Homogeneous Society", "Homogeneous Society",
    "Homogeneous Society", "Homogeneous Society",
    "Homogeneous Society", "Homogeneous Society",
    "Homogeneous Society", "Homogeneous Society",
    "Hypo-chaotic minority", "Hypo-chaotic minority",
    "Hypo-chaotic minority", "Hypo-chaotic minority",
    "Hypo-chaotic minority", "Hypo-chaotic minority",
    "Hypo-chaotic minority", "Hypo-chaotic minority",
    "Hypo-chaotic minority", "Hypo-chaotic minority",
    "Hypo-coupled minority", "Hypo-coupled minority",
    "Hypo-coupled minority", "Hypo-coupled minority",
    "Hypo-coupled minority", "Hypo-coupled minority",
    "Hypo-coupled minority", "Hypo-coupled minority",
    "Hypo-coupled minority", "Hypo-coupled minority"
  )
  # dist.a <- matrix(0, num.persons, nrow = num.persons)
  system.time(
    for(i in 1:num.persons){
      m1 <- l.now[[i]]$m %>% my_seriate()
      print(i)
      for(j in i:(num.persons)){
        m2 <- l.now[[j]]$m %>% my_seriate()
        dist.m[j,i] <- 1 - MatrixCorrelation::PSI(m1,m2)
        # dist.a[j,i] <- a %>% dist() %>% as.matrix()
        print(j)
      }
    }
  )
  
}


l.a <- list()
l.extracted <- list()
system.time(
for(b in brain_files_1e6){
  load(paste0(path.to.brains,"/",b))
  k <- brain_case %>% extract_now()
    k$owner %>% print()
  m <- k$a %>% 
    t() %>%  as.matrix() %>% 
    dist() %>% as.matrix()
  if(brain_case@now$activities %>% is.na() %>% sum()){
    m[,] <- 0
  }
  s <- m %>% seriate()
  l.a[[paste(k$verbal, k$owner, sep = "_")]] <- m[s[[1]], s[[2]]]
  l.extracted[[paste(k$verbal, k$owner, sep = "_")]] <- brain_case %>% extract_now()
  
}
)

o <- extract_now(b)

for (i in 1:25) {
  ind <- ((i - 1) * 20 + 1):(i * 20)
  paste(ind[[1]], "to", ind[[20]], "is being read now") %>% print()
  this.brain_location <- brain_locations[ind]
  
  Sys.time() %>% print()
  system.time(snp <- this.brain_location %>%
                ldply(extract_brains,
                      snapshots =  seq(5e3, 100e3, 5e3))) %>% print()
  
  save_vars(
    list.of.vars = "snp",
    prefix = paste0("snp_", ind[[1]], "-", ind[[20]]),
    path = path.to.snp
  )
}


dist.a <- matrix(0, num.persons, nrow = num.persons)
system.time(
  for(i in 1:50){
    print(i)
    for(j in i:(num.persons)){
      dist.a[j,i] <- 1 - MatrixCorrelation::PSI(l.a[[i]],l.a[[j]])
      print(j)
    }
  }
)


for(a in 1:length(l.a)){
  pimage(l.a[[a]], main = names(l.a)[a], key=FALSE)
}






# load("./data/snp.lean_all_5k_20190713_1356.RData")

a <- brain_case@now$activities %>% t()
brain_case@parameters$params.eps_a
d <- a %>% dist() %>% as.matrix()



d <- dist.a %>% as.matrix()
s <- d %>% seriate()
pimage(d,s)



c.d <- cluster::diana(d)
pltree(c.d, cex = 0.6, hang = -1,
       main = "diana activities")
rect.hclust(c.d, k = 5, border = 2:5)

c.a <- cluster::agnes(d)
pltree(c.a, cex = 0.6, hang = -1,
       main = "agnes activities")
rect.hclust(c.a, k = 5, border = 2:5)


fviz_dist(as.dist(d), gradient = list(low = "white", high = "black"))

# save_vars("l.extracted", )
rm(brain_case)
save_vars("l.extracted", prefix = "last-snapshot-conn-activity", path = "data-pc")
