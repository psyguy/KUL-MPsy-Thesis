# rm(list=ls())

load("data-pc/snp_only-1e+6_20190723_1404.RData")
source("functions/functions_my.R")

make.df <- function(input.df, col){
  
  d <- input.df %>% 
    select("Partition",
           "Owner",
           "Verbal.Description",
           col)
  
  output.df <- NULL
  for(p in 1:4){
    this.df <- select(d, -col)[p,]
    vec <- pull(d, col)[[p]]
    len <- length(vec)
    x <- this.df[rep(seq_len(nrow(this.df)),
                     each = len),]
    
    output.df <- output.df %>%
      rbind(
        cbind(x, vec, c(1:len))
      )
  }
  
  colnames(output.df) <- c("Partition",
                           "Owner",
                           "Verbal.Description",
                           col,
                           "Club Size")
  output.df %>% return()
}


rc_norm <- tmp$adj.mat.vect %>% 
  vec2mat() %>% 
  graph_from_adjacency_matrix() %>% 
  rich_club_norm(200)

snp.new <- snp %>% 
  filter(Partition == "whole") %>% 
  mutate(fn = paste(Verbal.Description, Owner, sep="_")) %>%
  arrange(fn) %>%
  mutate(fn = paste0(rep(c(3,2,1,4,5),each=10), "_",fn)) %>% 
  arrange(fn) %>% 
  cbind(family.codes = families, num =1:50) %>% 
  mutate(newcodes = paste0(num,"_",family.codes,rep(c(1:10), 5), "_", Owner)) %>% 
  filter(!(num %in%c(12,13,41,43))) # removing the dead brains
  




sig.level <- 0.05
significants <- data.frame(x = rc_norm$k[rc_norm$p < sig.level],
                           y = rc_norm$norm[rc_norm$p < sig.level])

  tmp <- snp %>%
    filter(Owner == name.this.owner) %>%
    filter(Rewiring == 1e+6) %>%
    mutate(`Vertex Betweenness` = map(adj.mat.vect, netmeas_bt,
                                      edge_betweenness = FALSE)) %>%
    mutate(`Edge Betweenness` = map(adj.mat.vect, netmeas_bt,
                                    edge_betweenness = TRUE)) %>%
    # mutate(`Centrality` = map(adj.mat.vect, netmeas_rc)) %>%
    mutate(`Rich Club` = map(adj.mat.vect, netmeas_rc))

m <- snp.new$adj.mat.vect %>% vec2mat() %>% graph_from_adjacency_matrix(mode="undirected")

rc.new <- function(x, N = 10){
  rco <- x %>%
    vec2mat() %>%
    graph_from_adjacency_matrix(mode="undirected") %>% 
    rich_club_norm(N=N)
  o <- cbind(rco$k,rco$norm, rco$p)
  
  colnames(o) <- c("Club Size", "Rich Club", "sig.level")
  return(o)
}


# rc <- NULL
# for(i in 1:46){
#   print(i)
#   rc.tmp <- snp.new$adj.mat.vect[[i]] %>% rc.new(200)
#   rc.name <- rep(snp.new$newcodes, nrow(rc.tmp))
#   rc <- rbind(rc, cbind(rc.name,rc.tmp))
# }
# save_vars(c("rc", "snp.new"), prefix = "NormalizedRC-onlyWhole")


load("data/NormalizedRC-onlyWhole_20200103_2046.RData")
rc <- rc %>%
  as.data.frame

rc$`Club Size` <- rc$`Club Size` %>% as.character() %>% as.numeric()
rc$`Rich Club` <- rc$`Rich Club` %>% as.character() %>% as.numeric()
rc$sig.level <- rc$sig.level %>% as.character() %>% as.numeric()

tmp <- rc %>% 
  mutate(n = gsub("_.*","", rc.name) %>% as.numeric) %>%
  mutate(family.num = (x-1)%/%10 + 1) %>% 
  filter(family.num == 1)
tmp$family.num <- tmp$family.num %>% as.factor()
tmp$n <- tmp$n %>% as.factor()


p <- ggplot(data = tmp,
            aes(x = `Club Size`,
                y = `Rich Club`,
                colour = n)) +
  geom_line(size = 1.5, alpha = 0.8) +
  # scale_colour_manual(values = c(colors$inter, colors$majo,
  #                                colors$mino, colors$whole)) +
  ggplot2::xlim(0,150)




Rich.Club.150 <- 
    tmp %>%
    ggplot(aes(x = `Club Size`,
               y = `Rich Club`, 
               colour = family.num)) +
    geom_line(size = 1.5, alpha = .8) +
    scale_colour_manual(values = c(colors$inter, colors$majo,
                                   colors$mino, colors$whole)) +
    theme(legend.position = "none") +
    geom_point(aes(x=x, y=y),
               data = significants,
               colour="black",
               size = 2.2,
               shape = 18)
  ggplot2::xlim(0, 150) 




x <- m %>% rich_club_norm(N=10)



# garbage -----------------------------------------------------------------

# rc_norm <- function(name.this.owner = NULL,
#                                  this.Verbal.Description = NULL,
#                                  this.Partition = NULL,
#                                  snp,
#                                  colors = list(
#                                    bg = "white",
#                                    mino = "deepskyblue3",
#                                    majo = "orangered",
#                                    inter = "olivedrab2",
#                                    whole = "dimgray"),
#                                  path.fig = "figures"){
#   
#   tmp <- snp %>%
#     filter(Owner == name.this.owner) %>% 
#     filter(Rewiring == 1e+6) %>% 
#     mutate(`Vertex Betweenness` = map(adj.mat.vect, netmeas_bt,
#                                       edge_betweenness = FALSE)) %>%
#     mutate(`Edge Betweenness` = map(adj.mat.vect, netmeas_bt,
#                                     edge_betweenness = TRUE)) %>%
#     # mutate(`Centrality` = map(adj.mat.vect, netmeas_rc)) %>%
#     mutate(`Rich Club` = map(adj.mat.vect, netmeas_rc))
#   
#   tmp$`Rich Club`[[2]] <- tmp$`Rich Club`[[2]]/tmp$`Rich Club`[[1]]
#   tmp$`Rich Club`[[3]] <- tmp$`Rich Club`[[3]]/tmp$`Rich Club`[[1]]
#   tmp$`Rich Club`[[4]] <- tmp$`Rich Club`[[4]]/tmp$`Rich Club`[[1]]
#   
#   rc_norm <- tmp$adj.mat.vect %>% 
#     vec2mat() %>% 
#     graph_from_adjacency_matrix() %>% 
#     rich_club_norm(200)
#   tmp$`Rich Club`[[1]] <- rc_norm$norm
#   
#   sig.level <- 0.05
#   significants <- data.frame(x = rc_norm$k[rc_norm$p < sig.level],
#                              y = rc_norm$norm[rc_norm$p < sig.level])
#   
#   Rich.Club.150 <- 
#     tmp %>%
#     make.df("Rich Club") %>%
#     # filter(Partition=="whole") %>% 
#     ggplot(aes(x = `Club Size`,
#                y = `Rich Club`, 
#                colour = Partition)) +
#     geom_line(size = 1.5, alpha = .8) +
#     scale_colour_manual(values = c(colors$inter, colors$majo,
#                                    colors$mino, colors$whole)) +
#     theme(legend.position = "none") +
#     geom_point(aes(x=x, y=y),
#                data = significants,
#                colour="black",
#                size = 2.2,
#                shape = 18)
#   ggplot2::xlim(0, 150) 
#   
#   
#   `Vertex Betweenness` <- tmp %>%
#     make.df("Vertex Betweenness") %>%
#     ggplot(aes(`Vertex Betweenness`)) +
#     geom_density(aes(fill = Partition), alpha = 0.6) +
#     scale_fill_manual(values = c(colors$inter, colors$majo,
#                                  colors$mino, colors$whole)) +
#     theme(legend.position = "none") +
#     ggplot2::xlim(0, 1000)
#   
#   
#   `Edge Betweenness` <- tmp %>%
#     make.df("Edge Betweenness") %>%
#     ggplot(aes(`Edge Betweenness`)) +
#     geom_density(aes(fill = Partition), alpha = 0.6) +
#     scale_fill_manual(values = c(colors$inter, colors$majo,
#                                  colors$mino, colors$whole)) +
#     theme(legend.position = "none") +
#     ggplot2::xlim(0, 150)
#   
#   
#   figure <- ggarrange(Rich.Club.150,
#                       `Vertex Betweenness`,
#                       `Edge Betweenness`,
#                       ncol = 1, nrow = 3,
#                       common.legend = TRUE,
#                       legend = "bottom"#ifelse(is.null(name.this.owner),
#                       # "bottom", "bottom")#"none", "bottom")
#   )
#   
#   vd <- snp %>%
#     filter(Owner == name.this.owner) %>%
#     pull(Verbal.Description) %>% 
#     as.character()
#   
#   title <- paste0("Final Rich Club and Betweenness of ",
#                   name.this.owner,
#                   " (", tolower(vd[[1]]), ")")
#   
#   
#   
#   pf <- path.fig
#   if(substr(pf, nchar(pf), nchar(pf))!="/") path.fig <- paste0(path.fig, "/")
#   file.name <- paste0(path.fig, title)
#   
#   annotate_figure(figure,
#                   top = text_grob(label =  paste0("",
#                                                   "Final Rich Club and Betweenness values of",
#                                                   "\n",
#                                                   name.this.owner,
#                                                   " (",
#                                                   tolower(vd[1]),
#                                                   ")"
#                   ),
#                   size = 25,
#                   family = "Times"
#                   )
#   ) %>%
#     graph2pdf(height = 15, width = 15,
#               file = paste0(file.name,".pdf")
#     )
#   
# }
# 
# netmeas_rc <- function(v, non.normalized = TRUE, k = c(1:150)){
#   m_ <- v %>% 
#     vec2mat()
#   g_ <- m_ %>% 
#     graph_from_adjacency_matrix()
#   
#   o1 <- 
#     sapply(k, function(x) brainGraph::rich_club_coeff(g_,
#                                                       k=x)$phi %>% as.numeric())
#   # if(sum(m_) != 10400 | non.normalized == TRUE) 
#   o1 %>% return()
#   # if(sum(m_) == 10400 & non.normalized == FALSE){
#   # z <-  g_ %>% rich_club_norm(200)
#   # o <- z$norm
#   # o %>% return()}
# }



