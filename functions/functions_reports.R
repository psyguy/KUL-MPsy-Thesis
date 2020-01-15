# rm(list = ls())

source("./functions/functions_partition.R")


# Plotting the connectivity matricies and networks at snapshots -----------


reports_netviz <- function(brain_case,
                           snapshots = c(50000, 100000),
                           width.column.report = 2240,
                           colors = list(
                             bg = "white",
                             mino = "deepskyblue3", #"blue2"
                             majo = "orangered", #"brown3"
                             inter = "olivedrab2", #"green4"
                             whole = "azure4"
                           )){
  
  
  for(rew.base in snapshots){
    
    rew <- rew.base + brain_case@age$rewires - 100001
    
    paste("Working on", brain_case@name, "at rewiring", rew) %>% print()
    m <- brain_case@history$mat.connectivity[[rew]]
    g <- m %>% graph_from_adjacency_matrix(mode="undirected")
    
    
    # rewrew <- ifelse(brain_case@age$rewires>100001,
    #                  rew + 100000,
    #                  rew)
    title <- brain_case@name %>% paste("after", rew, "rewirings")
    
    V(g)$partition <- c(rep("minority", 50),
                        rep("majority", 250))
    
    # select edges and set color 
    E(g)[V(g)[partition == "minority"] %--% V(g)[partition == "minority"]]$color <- colors$mino
    E(g)[V(g)[partition == "minority"] %--% V(g)[partition == "majority"]]$color <- colors$inter
    E(g)[V(g)[partition == "majority"] %--% V(g)[partition == "majority"]]$color <- colors$majo
    
    # plot
    set.seed(1)
    paste0(title, "_network", ".png") %>%
      png(width = width.column.report*2, height = width.column.report*2, res = 200)
    g %>% plot(vertex.size = 0,
               vertex.label = NA,
               edge.width = 1,
               edge.curved= 0.5
    )
    dev.off()
    
    m[1:50,1:50] <- m[1:50,1:50]*3
    m[1:50,51:300] <- m[1:50,51:300]*2
    m[51:300,1:50] <- m[51:300,1:50]*2
    
    width.column.report <- 2240
    col.pimage <- c(colors$bg,
                    colors$majo,
                    colors$inter,
                    colors$mino)
    
    
    paste0(title, "_unserialized", ".png") %>%
      png(width = width.column.report, height = width.column.report, res = 400)
    pimage(m,
           col = col.pimage,
           key = FALSE
    )
    dev.off()
    
    
    paste0(title, "_serialized", ".png") %>%
      png(width = width.column.report, height = width.column.report, res = 400)
    pimage(m,
           seriate(m),
           col = col.pimage,
           key = FALSE
    )
    dev.off()
    
    # putting adj mats next to each other
    gl <- rep(title,2) %>%
      paste0(c("_unserialized.png", "_serialized.png")) %>% 
      as.list() %>%
      lapply(png::readPNG) %>% 
      lapply(grid::rasterGrob)
    
    paste0(title, "_connectivities", ".png") %>%
      png(width = width.column.report*2, height = width.column.report, res = 400)
    gr.con <- gridExtra::grid.arrange(grobs=gl, ncol=2,
                                      padding = unit(1, "point")
    )
    dev.off()
    
    
    conns_and_net <- rep(title,2) %>%
      paste0(c("_connectivities.png", "_network.png")) %>% 
      as.list() %>%
      lapply(png::readPNG) %>% 
      lapply(grid::rasterGrob)
    require(grid)
    paste0(title, "", ".png") %>%
      png(width = width.column.report*2, height = width.column.report*3, res = 400)
    gr.net <- gridExtra::grid.arrange(grobs=conns_and_net, nrow=2,
                                      top = textGrob(paste("\n",
                                                           brain_case@name,
                                                           "after",
                                                           rew,
                                                           "rewirings"),
                                                     gp=gpar(fontsize=30,font=8)
                                      )
    )
    dev.off()
  }
}


reports_coefs.extract <- function(pattern,
                                  sampled.path = "./data/"){
  sampled.names <- list.files(path = sampled.path, pattern = "*.RData")
  r.this <- sampled.names[grepl(pattern, sampled.names)]
  r.this <- r.this[grepl("_g-", r.this)]
  r.this <- r.this[!grepl("coefs.wb.", r.this)]
  
  coefs.wb <- NULL
  t <- Sys.time()
  for(sampled in r.this){
    # rm(brain_case)
    load(paste(sampled.path, sampled, sep = ""))
    brain_case@name %>% print()
      for(i in 1:(brain_case@age$rewires-1)){
        m <- brain_case@history$mat.connectivity[[i]]
        if(is.null(m)) next
        if(!(i%%5000)) paste("Adding coefs of rewiring", i,
                             "for", brain_case@name) %>% print()
        coefs.wb <- coefs.wb %>% rbind(m %>% netmeas_wbcoefs(parameters = brain_case@parameters,
                                                                 name = brain_case@name,
                                                                 rewiring = i))
      }
  }
  (Sys.time()-t) %>% paste("for", pattern) %>% print()

  save_vars("coefs.wb", prefix = paste0("coefs.wb.hpc_", pattern))
}
