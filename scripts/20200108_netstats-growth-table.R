
source("./functions/functions_partition.R")

# making random graphs and calculating mean values ------------------------

r.seed <- 1:100
families <<- c("OC", "OT", "BL", "UT", "UC") %>% rep(each = 10)

ran.graphs <- r.seed %>% 
  map(~make_random_graph(size = 300,
                         num.links = 5200,
                         seed = .))

quick.coefs <- function(m_){
  g_ <- m_ %>% graph_from_adjacency_matrix(mode = "undirected")
  
  c_ <- m_ %>% my_clustceof()
  e_ <- m_ %>% netmeas_efficiency()
  modu_ <- g_ %>% cluster_fast_greedy() %>% modularity()
  pl_ <- g_ %>% average.path.length(unconnected = TRUE)
  
  
  s_ <- c_*e_
  
  o <- c(Clustering = c_,#/c_0
         Small.World = s_,
         Modularity = modu_,
         Average.Path.Length = pl_,
         Assortativity = g_ %>% assortativity.degree() %>% as.numeric()
  ) %>% return()
  
}

r.g_whole <- ran.graphs
r.g_mino <- ran.graphs %>% map(~.[1:50,1:50])
r.g_majo <- ran.graphs %>% map(~.[51:300,51:300])
r.g_inter <- r.g_whole %>%
  map(function(x){
    x[1:50,1:50] <- 0
    x[51:300,51:300] <- 0
    return(x)
  })
system.time(rand.means <- ran.graphs )


r.m_whole <- r.g_whole %>% ldply(quick.coefs) %>% colMeans() %>% abs()
r.m_mino <- r.g_mino %>% ldply(quick.coefs) %>% colMeans() %>% abs()
r.m_majo <- r.g_majo %>% ldply(quick.coefs) %>% colMeans() %>% abs()
r.m_inter <- r.g_inter %>% ldply(quick.coefs) %>% colMeans() %>% abs()

rand.means <- rbind(r.m_whole,r.m_mino,r.m_majo,r.m_inter) %>% as.data.frame()

load("./data-pc/snp.lean_all_5k_20190713_1356.RData")



snp.new <- snp.lean

snp.new$Clustering <- snp.new$Clustering/rep(rand.means$Clustering,10000)
snp.new$Small.World <- snp.new$Small.World/rep(rand.means$Small.World,10000)
snp.new$Modularity <- snp.new$Modularity/rep(rand.means$Modularity,10000)
snp.new$Average.Path.Length <- snp.new$Average.Path.Length/rep(rand.means$Average.Path.Length,10000)
# snp.new$Assortativity <- snp.new$Assortativity/rep(rand.means$Assortativity,10000)
snp.new$Edge.Density <- snp.new$Edge.Density/rep(0.2329177,40000)




snp.new <- snp.new %>% 
  # filter(Partition == "whole") %>% 
  mutate(fn = paste(Verbal.Description, Owner, sep="_")) %>%
  arrange(fn) %>%
  mutate(fn = paste0(rep(c(3,2,1,4,5),each=10), "_",fn)) %>% 
  arrange(fn) %>% 
  cbind(family.codes = families, num =1:50) %>% 
  mutate(newcodes = paste0(num,"_",family.codes,rep(c(1:10), 5), "_", Owner)) %>% 
  mutate(famcode.num = paste0(family.codes, rep(c(1:10), 5)))

# snp.new$family.codes <- snp.new$family.codes %>% as.factor()

netstats.growth <- snp.new %>% 
  filter(!(famcode.num %in% c("OT2", "OT3", "UC1", "UC3"))) %>% 
  filter(Rewiring>59999) %>%
  group_by(family.codes,Partition) %>% 
  dplyr::summarize_at(c("Clustering", "Average.Path.Length",
                        "Small.World", "Modularity",
                        "Assortativity", "Edge.Density"),
                      ~paste0(round(mean(.),digits = 2),
                              "(",
                              round(sd(.),digits = 2),
                              ")")
                      )

colnames(netstats.growth) <- c("Condition",
                               "Partition",
                               "Clustering coefficient",
                               "Average path length",
                               "Small-worldness",
                               "Modularity",
                               "Assortativity",
                               "Edge density") %>% 
  as.data.frame() #%>% select(-Partition)

save(netstats.growth, file = "data/netstats.growth.RData")

