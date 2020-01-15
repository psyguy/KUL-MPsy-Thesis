# rm(list = ls())

source("./functions/functions_partition.R")


# Plotting the connectivity matricies and networks at snapshots -----------

extract_partitions <- function(m){
  minority <- m[1:50,1:50]
  majority <- m[51:300,51:300]
  # removing the within group connections
  interpartition <- m
  interpartition[1:50,1:50] <- interpartition[51:300,51:300] <- 0
  m.list <- list(whole = m,
                 minority = minority,
                 majority = majority,
                 interpartition = interpartition)
  m.list %>% return()
}



extract_identifiers <- function(b){
  
  ## debug
  # b <- brain_case
  
  Owner <- b@name #%>% rep(l_)
  Seed <- b@parameters$seed #%>% rep(l_)
  Round <- b@parameters$round
  p_d <- b@parameters$params.eps_a
  Epsilon.Proportion <- paste0("(",p_d[1],", ",p_d[2],", ",p_d[3],")")
  a.Proportion <- paste0("(",p_d[4],", ",p_d[5],", ",p_d[6],")")
  
  vd <- "Homogeneous Society"
  if(Epsilon.Proportion == "(0, 5, 1)") vd <- "Hyper-coupled minority"
  if(Epsilon.Proportion == "(1, 5, 0)") vd <- "Hypo-coupled minority"
  if(a.Proportion == "(0, 5, 1)") vd <- "Hyper-chaotic minority"
  if(a.Proportion == "(1, 5, 0)") vd <- "Hypo-chaotic minority"
  Verbal.Description <- vd
  
  identifiers <- data.frame(Owner,
                            Seed,
                            Round,
                            Rewiring = 0,
                            Verbal.Description,
                            Epsilon.Proportion,
                            a.Proportion)
  identifiers <- lapply(1:4, function(x) identifiers)
  names(identifiers) <- c("whole",
                          "minority",
                          "majority",
                          "interpartition")
  identifiers  %>% return()
}



extract_coefs <- function(m,
                          identifiers = NULL){
  
  g_ <- m %>% graph_from_adjacency_matrix(mode = "undirected")
  
  Clustering <- m %>% my_clustceof()
  Efficiency <- m %>% netmeas_efficiency()
  Small.World <- Clustering*Efficiency
  Modularity <- g_ %>% cluster_fast_greedy() %>% modularity()
  Assortativity <- g_ %>% assortativity.degree() %>% as.numeric()
  Rich.Club <- brainGraph::rich_club_coeff(g_)$phi %>% as.numeric()
  Average.Path.Length <- g_ %>% average.path.length(unconnected = TRUE)
  Edge.Density <- m %>% sum() %>% sum()
  
  coefs <- data.frame(Clustering, 
                      Small.World,
                      Modularity,
                      Assortativity,
                      Rich.Club,
                      Average.Path.Length,
                      Edge.Density,
                      Efficiency
  )
  
  if(is.data.frame(identifiers)) coefs <- cbind(identifiers[1,],
                                                coefs)
  
  coefs %>% return()
}



extract_id.n.mats <- function(m_raw, identifiers){
  
  m_partitions <- m_raw %>%
    extract_partitions()
  coefs_partitioned <- m_partitions %>%
    map(extract_coefs)
  
  output <- identifiers %>%
    map2(coefs_partitioned, cbind) %>% 
    bind_rows(.id = "Partition") %>% 
    mutate(adj.mat.vect = map(m_partitions, mat2vec))
  
  # normalizing the edge densities
  denom <- c(300, 50, 250)
  denom <- (denom*(denom-1)/2) %>% c(50*250)
  output$Edge.Density <- output$Edge.Density / rep(denom, nrow(output)/4)
  output %>% return()
}


extract_brains <- function(b_loc,
                           snapshots = c(50e3, 100e3),
                           incl.activities = FALSE){
  ## debug
  # b_loc <- this.brain_location
  # snapshots <- 100e3
  
  library(purrr)
  
  load(b_loc)
  age_ <- brain_case@age$rewires - 1
  
  brain_case@name %>% paste("being extracted at age", age_) %>% print()
  
  l.m <- brain_case@history$mat.connectivity
  l.m[1] <- NULL
  l.m.compact <- l.m %>% compact()
  
  identifiers <- brain_case %>% extract_identifiers()
  
  output <- l.m.compact[snapshots/200] %>%
    ldply(extract_id.n.mats, identifiers)
  
  output$Rewiring <- rep((age_ - 100e3 + snapshots), each = 4)
  
  output %>% return()
}


extract_plotnet <- function(m,
                            title = "Network",
                            colors = list(
                              bg = "white",
                              mino = "deepskyblue3", #"blue2"
                              majo = "orangered", #"brown3"
                              inter = "olivedrab2", #"green4"
                              whole = "azure4"),
                            save = TRUE,
                            first.add = FALSE,
                            curve = 0,
                            ps = 3,
                            path.fig = "figures"
                            ){
  if(is.vector(m)) m <- m %>% vec2mat()
  pf <- path.fig
  if(substr(pf, nchar(pf), nchar(pf))!="/") path.fig <- paste0(path.fig, "/")
  
  g <- m %>% graph_from_adjacency_matrix(mode="undirected")
  
  V(g)$partition <- c(rep("minority", 50),
                      rep("majority", 250))
  
  
  #Set the margin size (small margins)
  par(mar = rep(0.05, 4))
  vertex.color <- c("skyblue", "salmon")[1 + (V(g)$partition == "majority")]
  vertex.size <- 2
  ## select edges and set color and plot one by one
  # First remove them, then add majority, plot majority 
  E(g)$color <- NA
  E(g)[V(g)[partition == "majority"] %--% V(g)[partition == "majority"]]$color <- colors$majo
  set.seed(1)
  g %>% plot(vertex.size = vertex.size*sqrt(ps)/2,
             vertex.color = vertex.color,
             add = first.add,
             vertex.label = NA,
             edge.width = 1*ps,
             edge.curved= curve*ps)
  
  E(g)$color <- NA
  E(g)[V(g)[partition == "minority"] %--% V(g)[partition == "majority"]]$color <- colors$inter
  set.seed(1)
  g %>% plot(vertex.size = vertex.size*sqrt(ps)/2,
             vertex.color = vertex.color,
             add=TRUE,
             vertex.label = NA,
             edge.width = 1*ps,
             edge.curved= curve*ps)
  
  
  E(g)$color <- NA
  E(g)[V(g)[partition == "minority"] %--% V(g)[partition == "minority"]]$color <- colors$mino
  set.seed(1)
  g %>% plot(vertex.size = vertex.size*sqrt(ps)/2,
             vertex.color = vertex.color,
             add=TRUE,
             vertex.label = NA,
             edge.width = 1*ps,
             edge.curved= curve*ps)
  
  if(save) graph2pdf(height = 10*ps, width = 10*ps,
                     file = paste0(path.fig, title, "_network"))
}


extract_plotcon <- function(m,
                            title = "Connectivity",
                            colors = list(
                              bg = "white",
                              mino = "deepskyblue3", #"blue2"
                              majo = "orangered", #"brown3"
                              inter = "olivedrab2", #"green4"
                              whole = "dimgray"),
                            save = TRUE,
                            ps = 3,
                            path.fig = "figures"){
  if(is.vector(m)) m <- m %>% vec2mat
  pf <- path.fig
  if(substr(pf, nchar(pf), nchar(pf))!="/") path.fig <- paste0(path.fig, "/")
  file.name <- paste0(path.fig, title)
  
  s <- m %>% seriate(method="PCA_angle")
  m[1:50,1:50] <- m[1:50,1:50]*3
  m[1:50,51:300] <- m[1:50,51:300]*2
  m[51:300,1:50] <- m[51:300,1:50]*2
  
  col.pimage <- c(colors$bg, colors$majo,
                  colors$inter, colors$mino)
  
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow = 1, ncol = 2)))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  
  pimage(m,
         col = col.pimage,
         key = FALSE,
         newpage = FALSE)
  
  # if(save) graph2pdf(height = 5, width = 5,
  #                    file = paste0(path.fig, title, "_unserialized"))
  
  upViewport(1)
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  
  pimage(m,
         s,
         col = col.pimage,
         key = FALSE,
         newpage = FALSE)
  
  upViewport(1)
  popViewport(0)
  
  if(save) graph2pdf(height = 5*ps, width = 5*2*ps,
                     file = paste0(path.fig, title, "_connectivities"))
  
}


extract_plotglue <- function(title = "Someone",
                             path.fig = "figures"){
  pf <- path.fig
  if(substr(pf, nchar(pf), nchar(pf))!="/") path.fig <- paste0(path.fig, "/")
  file.name <- paste0(path.fig, title)
  
  panel.upper.l <- file.name %>% paste0("_unserialized.pdf") %>% image_read_pdf()
  panel.upper.r <- file.name %>% paste0("_serialized.pdf") %>% image_read_pdf()
  panel.lower <- file.name %>% paste0("_network.pdf") %>% image_read_pdf()
  
  whole <- c(panel.upper.l, panel.upper.r) %>% 
    image_append() %>% 
    c(panel.lower) %>% 
    image_append(stack = TRUE)
  
  image_write(whole,"a.pdf")
  
  graph2pdf(whole,
            height = 15, width = 10,
            file = paste0("Profile of", title))
  
}


extract_plotcoefs.single <- function(chosen.coef,
                                     name.this.owner = NULL,
                                     this.Verbal.Description = NULL,
                                     this.Partition = NULL,
                                     snp,
                                     colors = list(
                                       bg = "white",
                                       mino = "deepskyblue3",
                                       majo = "orangered",
                                       inter = "olivedrab2",
                                       whole = "dimgray")){
  if(!is.null(name.this.owner)){
    p <- ggplot(data = snp %>% filter(Owner == name.this.owner),
                aes(x = Rewiring,
                    y = !!ensym(chosen.coef),
                    colour = Partition)) +
      geom_line(size = 1.5, alpha = 0.8) +
      scale_colour_manual(values = c(colors$inter, colors$majo,
                                     colors$mino, colors$whole)) +
      ggplot2::ylim(min(0, min(snp[chosen.coef])),
                    max(snp[chosen.coef]))
    p %>% return()
  }
  
  p <- ggplot(data = snp %>%
                filter(Verbal.Description == this.Verbal.Description) %>%
                filter(Partition == this.Partition),
              aes(x = Rewiring,
                  y = !!ensym(chosen.coef),
                  colour = Owner)) +
    geom_line(size = 0.5, alpha = 0.8) +
    # scale_colour_manual(values = c(colors$inter, colors$majo,
    #                                colors$mino, colors$whole)) +
    theme(legend.position = "none")+
    ggplot2::ylim(min(0, min(snp[chosen.coef])),
                  max(snp[chosen.coef]))
  p %>% return()
  
}


extract_plotcoefs.glued <- function(name.this.owner = NULL,
                                    this.Verbal.Description = NULL,
                                    this.Partition = NULL,
                                    snp,
                                    coef.range = c(9:16), # which coefficients to plot
                                    path.fig = "figures"){
  coef.names <- names(snp)[coef.range] %>% as.list()
  list.of.plots <- coef.names %>% 
    llply(extract_plotcoefs.single,
          name.this.owner = name.this.owner,
          this.Verbal.Description = this.Verbal.Description,
          this.Partition = this.Partition,
          snp = snp)
  
  figure <- ggarrange(plotlist = list.of.plots,
                      ncol = 2, nrow = 3,
                      common.legend = TRUE,
                      legend = ifelse(is.null(name.this.owner),
                                      "bottom", "bottom")#"none", "bottom")
                      )
  
  if(is.null(name.this.owner)){
    vd <- snp %>%
      filter(Verbal.Description == this.Verbal.Description) %>%
      filter(Partition == this.Partition) %>% 
      pull(Verbal.Description) %>% 
      as.character()
  }else{
    vd <- snp %>%
      filter(Owner == name.this.owner) %>%
      pull(Verbal.Description) %>% 
      as.character()
  }
  
  this.Partition <- ifelse(is.null(this.Partition),
                           "",
                           paste0(this.Partition, " "))
  
  title <- paste0("Network statistics of ",
                  name.this.owner,
                  this.Partition,
                  "(", tolower(vd[[1]]), ")")
  
  pf <- path.fig
  if(substr(pf, nchar(pf), nchar(pf))!="/") path.fig <- paste0(path.fig, "/")
  file.name <- paste0(path.fig, title)
  
  
  annotate_figure(figure,
                  top = text_grob(label =  paste0("",
                                                  "Network statistics of ",
                                                  name.this.owner,
                                                  this.Partition,
                                                  "\n (",
                                                  tolower(vd[1]),
                                                  ")"),
                                  size = 25, family = "Times"
                                  )
                  ) %>%
    graph2pdf(height = 15, width = 15,
              file = file.name
              )
  
}


extract_plot_rc.btwn <- function(name.this.owner = NULL,
                                 this.Verbal.Description = NULL,
                                 this.Partition = NULL,
                                 snp,
                                 colors = list(
                                   bg = "white",
                                   mino = "deepskyblue3",
                                   majo = "orangered",
                                   inter = "olivedrab2",
                                   whole = "dimgray"),
                                 path.fig = "figures"){
  
  tmp <- snp %>%
    filter(Owner == name.this.owner) %>% 
    filter(Rewiring == 1e+6) %>% 
    mutate(`Vertex Betweenness` = map(adj.mat.vect, netmeas_bt,
                                      edge_betweenness = FALSE)) %>%
    mutate(`Edge Betweenness` = map(adj.mat.vect, netmeas_bt,
                                    edge_betweenness = TRUE)) %>%
    # mutate(`Centrality` = map(adj.mat.vect, netmeas_rc)) %>%
    mutate(`Rich Club` = map(adj.mat.vect, netmeas_rc))
  
  tmp$`Rich Club`[[2]] <- tmp$`Rich Club`[[2]]/tmp$`Rich Club`[[1]]
  tmp$`Rich Club`[[3]] <- tmp$`Rich Club`[[3]]/tmp$`Rich Club`[[1]]
  tmp$`Rich Club`[[4]] <- tmp$`Rich Club`[[4]]/tmp$`Rich Club`[[1]]
  
  rc_norm <- tmp$adj.mat.vect %>% 
    vec2mat() %>% 
    graph_from_adjacency_matrix() %>% 
    rich_club_norm(200)
  tmp$`Rich Club`[[1]] <- rc_norm$norm
  
  sig.level <- 0.05
  significants <- data.frame(x = rc_norm$k[rc_norm$p < sig.level],
                             y = rc_norm$norm[rc_norm$p < sig.level])
    
  Rich.Club.150 <- 
    tmp %>%
    make.df("Rich Club") %>%
    # filter(Partition=="whole") %>% 
    ggplot(aes(x = `Club Size`,
               y = `Rich Club`, 
               colour = Partition)) +
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
  
  
  `Vertex Betweenness` <- tmp %>%
    make.df("Vertex Betweenness") %>%
    ggplot(aes(`Vertex Betweenness`)) +
    geom_density(aes(fill = Partition), alpha = 0.6) +
    scale_fill_manual(values = c(colors$inter, colors$majo,
                                 colors$mino, colors$whole)) +
    theme(legend.position = "none") +
    ggplot2::xlim(0, 1000)
  
  
  `Edge Betweenness` <- tmp %>%
    make.df("Edge Betweenness") %>%
    ggplot(aes(`Edge Betweenness`)) +
    geom_density(aes(fill = Partition), alpha = 0.6) +
    scale_fill_manual(values = c(colors$inter, colors$majo,
                                 colors$mino, colors$whole)) +
    theme(legend.position = "none") +
    ggplot2::xlim(0, 150)
  
  
  figure <- ggarrange(Rich.Club.150,
                      `Vertex Betweenness`,
                      `Edge Betweenness`,
                      ncol = 1, nrow = 3,
                      common.legend = TRUE,
                      legend = "bottom"#ifelse(is.null(name.this.owner),
                      # "bottom", "bottom")#"none", "bottom")
  )
  
    vd <- snp %>%
    filter(Owner == name.this.owner) %>%
    pull(Verbal.Description) %>% 
    as.character()
  
  title <- paste0("Final Rich Club and Betweenness of ",
                  name.this.owner,
                  " (", tolower(vd[[1]]), ")")
  

  
  pf <- path.fig
  if(substr(pf, nchar(pf), nchar(pf))!="/") path.fig <- paste0(path.fig, "/")
  file.name <- paste0(path.fig, title)
  
  annotate_figure(figure,
                  top = text_grob(label =  paste0("",
                                                  "Final Rich Club and Betweenness values of",
                                                  "\n",
                                                  name.this.owner,
                                                  " (",
                                                  tolower(vd[1]),
                                                  ")"
                  ),
                  size = 25,
                  family = "Times"
                  )
  ) %>%
    graph2pdf(height = 15, width = 15,
              file = paste0(file.name,".pdf")
    )
  
}


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
