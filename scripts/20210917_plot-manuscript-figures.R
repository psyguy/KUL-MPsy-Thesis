# Notes -------------------------------------------------------------------
##########################################################
######## plotting adj+net with PCA_angle seriation #######
######## 2021-09-17 ######################################
##########################################################
## 
## From commit 2aa903a on, the following snp data files should be used:
##    - snp_all_5k_20210628_1238.RData (zipped into three parts, for GitHub)
##    - snp.lean_all_5k_20210628_1241.RData
##    - snp_only-1e+6_20210628_1241.RData
## 
## Here I change names of SD/HD to LC/MC, which takes place at labeling
##########################################################


# Loading functions and snp -----------------------------------------------

library(librarian)
shelf(cowplot)
shelf(tikzDevice)
source("functions/functions_my.R")
source("functions/functions_netmeas.R")
source("./functions/functions_extract.R")
load("./data/snp_only-1e+6_20210628_1241.RData")
load("./data/snp.lean_all_5k_20210628_1241.RData")
load("./data/snp.rc_only-1e+6_20210701_2257.RData")


make.df <- function(input.df, variable){
  
  d <- input.df %>% 
    select("Partition",
           "ID",
           "Verbal.Description",
           all_of(variable))
  
  output.df <- NULL
  for(p in 1:4){
    this.df <- select(d, -variable)[p,]
    vec <- pull(d, variable)[[p]]
    len <- length(vec)
    x <- this.df[rep(seq_len(nrow(this.df)),
                     each = len),]
    
    output.df <- output.df %>%
      rbind(
        cbind(x, vec, c(1:len))
      )
  }
  
  colnames(output.df) <- c("Partition",
                           "ID",
                           "Verbal.Description",
                           variable,
                           "Club Size")
  output.df %>% return()
}


plot_labels <- function(label = "Whole",
                        size = 5,
                        angle = 0,
                        family = "CMU Serif"){
  ggplot() +
    annotate(geom = 'text',
             x=1, y=1,
             label = label,
             size = size,
             family = family,
             angle = angle) +
    theme_void()
}

# Change font in plots ----------------------------------------------------
# https://stackoverflow.com/a/51906008/2275986
shelf(showtext)

# Check the current search path for fonts
font_paths()    
#> [1] "C:\\Windows\\Fonts"

# List available font files in the search path
f <- font_files()    

# syntax: font_add(family = "<family_name>", regular = "/path/to/font/file")
font_add("CMU Classical Serif", "cmunci.ttf")
font_add("CMU Serif Upright Italic", "cmunui.ttf")
font_add("CMU Serif", "cmunrm.ttf")

font_families()

## automatically use showtext for new devices
showtext_auto() 


# Getting only whole partition --------------------------------------------

snp.whole <- snp %>% 
  filter(Partition == "whole")

all.owners <- snp.whole$fn %>% 
  unique() %>% 
  sort() %>% 
  as.character() %>%
  as.list()


# Saving connectivity and network plots (fig 1 & 4) -----------------------

# It saves plots separately (connectivity+seriated & networks) in separate PDFs
# Has to be combined manually
for(i in 1:nrow(snp.whole)){
  m <- snp.whole$adj.mat.vect[[i]] %>% vec2mat()
  title <- snp.whole$ID[[i]]
  print(paste (i, title))
  
  extract_plotcon(m, title)
  extract_plotnet(m, title)
  
}


# Plotting evolutions -----------------------------------------------------

plot_evolution <- function(df = snp.lean,
                           variable = "Modularity",
                           fam.code = "BL",
                           Partition_ = "whole",
                           scale.down = 1){
  # df <- snp.lean
  df <- df %>% 
    select(Partition,
           ID,
           fn,
           Rewiring,
           all_of(variable)
           ) %>% 
    filter(Partition == Partition_) %>%
    filter(grepl(fam.code, ID, fixed = TRUE))
  colnames(df)[ncol(df)] <- "value"
  
  df$Rewiring <- df$Rewiring/1000
  
  max.y <- 1
  
  if(variable == "Average.Path.Length") max.y <- 6
  if(variable == "Edge.Density") max.y <- 2
  
  varname <- gsub("\\.", " ", variable)
  if(variable == "Small.World") varname <- "Small-world Index"
  if(variable == "Clustering") varname <- "Clustering Coefficient"
  
  font.cmu.serif.upright.italic <- "CMU Serif Upright Italic"
  font.cmu.classical.serif <- "CMU Classical Serif"
  
  ev.plot <- df %>% 
    ggplot(
      aes(x = Rewiring,
        y = value,
        colour = ID)) +
    geom_line(size = 0.2*scale.down, alpha = .8) +
    # labs(x = "Rewiring",
    #      y = "Value"
    # ) +
    # theme_half_open() +
    theme_bw() +
    theme(text = element_text(size = 10, family = "CMU Serif"
                              ),
          # plot.title = element_text(hjust = 0.5,
          #                           size = ifelse(variable == "Edge.Density",
          #                                         20,
          #                                         0.01)),
          # plot.subtitle = element_text(size = ifelse(variable == "Edge.Density",
          #                                            5,
          #                                            0.01)),
          # axis.title.x = element_blank(),
          # axis.title.y = element_blank(),
          legend.position = "none"
    ) +
    labs(
         # title = ifelse(variable == "Edge.Density",
         #                stringr::str_to_title(Partition_),
         #                ""),
         # subtitle = "",
         x = ifelse(variable == "Assortativity", "Thousand Rewirings",""),
         y = ifelse(Partition_ == "whole", varname,"")
         ) +
    # annotate('text',
    #          x = 0,
    #          y = max.y,
    #          hjust = 0,
    #          vjust = 1,
    #          size = 3.5*scale.down,
    #          # family = font.cmu.serif.upright.italic,
    #          label = gsub("\\.", " ", variable)
    # ) +
  scale_x_continuous(breaks = seq(0, 1000, 250),
                     limits = c(0, 1000)) +
    scale_y_continuous(breaks = seq(0, max.y, max.y*0.2),
                       limits = c(0, max.y))
  return(ev.plot)
  
}

plot_evolution.all <- function(df = snp.lean,
                   fam.code = "BL",
                   Partition_ = "whole",
                   scale.down = .85){
  plot_grid(
    plot_labels(paste0("  ",
                      stringr::str_to_title(Partition_)
                      ),
                size = 8),
    plot_labels(""),
    plot_evolution(df,"Edge.Density",
                   fam.code, Partition_, scale.down),
    plot_evolution(df,"Clustering",
                   fam.code, Partition_, scale.down),
    plot_evolution(df,"Average.Path.Length",
                   fam.code, Partition_, scale.down),
    plot_evolution(df,"Small.World",
                   fam.code, Partition_, scale.down),
    plot_evolution(df,"Modularity",
                   fam.code, Partition_, scale.down),
    plot_evolution(df,"Assortativity",
                   fam.code, Partition_, scale.down),
    ncol = 1,
    rel_heights = c(1.2,0.5,6,6,6,6,6,6),
    byrow = TRUE)
}

spacer.margin <- 0.2
spacer <- 0.1

plot_ev.perfamily <- function(fam.code = "HC",
                              # col.title.height = 0.2,
                              # spacer.margin = spacer.margin,
                              # spacer = spacer,
                              df = snp.lean,
                              scale.down = .85){
  
  
  plot.ev.all <- plot_grid(
    plot_evolution.all(fam.code = fam.code, Partition_ = "whole"),
    plot_evolution.all(fam.code = fam.code, Partition_ = "majority"),
    plot_evolution.all(fam.code = fam.code, Partition_ = "minority"),
    ncol = 3,
    rel_widths = c(1,1,1))
  
  
  if(fam.code == "SD") fam.code <- "LC"
  if(fam.code == "HD") fam.code <- "MC"
  
  save_plot(paste0("evolution-",fam.code,".pdf"),
            plot.ev.all,
            base_height = 11.69,
            base_width = 8.27)
  
}

plot_ev.perfamily("BL")
plot_ev.perfamily("SD")
plot_ev.perfamily("HD")
plot_ev.perfamily("SC")
plot_ev.perfamily("HC")

tikz(file = "plot_test.tex", width = 8.27, height = 11.69)

print(plot.final)

dev.off()


# Calculating rich club coefficient ---------------------------------------

# removing the dead brains
snp.new <- snp %>%
  # filter(Partition != )
  filter(!(ID %in% c("HD2", "HD3", "SC1", "SC3")))


snp.rc <- NULL
system.time(
for(r in 1:nrow(snp.new)){
  print(paste("Now calculating",
              r,
              snp.new[r,]$ID,
              snp.new[r,]$Partition))
  system.time(
  d.tmp <- netmeas_richclub(snp.new[r,],
                            k.max = 200,
                            N.rand = 200)
  ) %>% print()
  if(is.null(snp.rc)) snp.rc <- d.tmp
  else snp.rc <- rbind(snp.rc,d.tmp)
}
)


save_vars("snp.rc", prefix = "snp.rc_only-1e+6")


# plotting RC values and significance levels ------------------------------

rc.res <- snp.rc

# rc.res$num <- rc.res$num %>% as.character() %>% as.factor()

plot_richclub <- function(rc.res = snp.rc,
                          fam.code = "BL",
                          Partition_ = "whole"){
  scale.down <- 0.5
  if(Partition_ == "whole") max.clubsize <- 85
  if(Partition_ == "majority") max.clubsize <- 65
  if(Partition_ == "minority") max.clubsize <- 40
  
# fam.code <- "SC"
# Partition_ <- "whole"

rc.this <- rc.res %>%
  filter(Partition == Partition_) %>%
  mutate(rc.p = `Normalized Rich Club`) %>% 
  filter(grepl(fam.code, Verbal.Description, fixed = TRUE))

# rc.this$`Normalized Rich Club`[rc.this$`Normalized Rich Club`==0] <- NA
rc.this$rc.p[is.nan(rc.this$p.value)] <- NA
rc.this$`Normalized Rich Club`[is.nan(rc.this$`Normalized Rich Club`)] <- NA
rc.this <- rc.this %>% na.omit()
rc.this$rc.p[rc.this$p.value>0.01] <- NA

font.cmu.serif.upright.italic <- "CMU Serif Upright Italic"
font.cmu.classical.serif <- "CMU Classical Serif"

rc.plot <- rc.this %>%  ggplot(
  aes(x = `Club Size`,
      y = `Normalized Rich Club`,
      colour = ID)) +
  geom_line(size = 1*scale.down, alpha = 0.45) +
  labs(x = ifelse(fam.code == "SC", "Club size",""),
       # y = ifelse(Partition_ == "whole", "Value","")
       y = ""
       # title = paste("Normalized Rich-Club for",
       #               fam.code,
       #               paste0("(", Partition_, ")"))
  ) +
  # theme_half_open() +
  theme_bw() +
  # labs(x = ifelse(fam.code == "SD", "Club size",element_blank()),
  #      y = ifelse(Partition_ == "whole", "Value",element_blank()),
  #      # title = paste("Normalized Rich-Club for",
  #      #               fam.code,
  #      #               paste0("(", Partition_, ")"))
  #      ) +
  theme(text = element_text(size = 12, family = font.cmu.classical.serif),
        # axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10,
                                  family = font.cmu.serif.upright.italic,
                                  hjust = 0.5,
                                  vjust = -12),
        # plot.subtitle = element_text(size = 15,
        #                              family = font.cmu.serif.upright.italic,
        #                              hjust = 0.01,
        #                              vjust = -18),
        legend.position = "none"
        ) +
    geom_point(aes(x = `Club Size`,
                 y = rc.p,
                 color = ID),
             data = rc.this,
             size = 0.5*scale.down,
             shape = 16) +
  geom_hline(yintercept = 1,
             size = 0.1,
             linetype = "dashed") +
    # annotate('text',
    #          x = 0,
    #          y = 4.2,
    #          hjust = 0,
    #          vjust = 1,
    #          size = 10,
    #          family = font.cmu.serif.upright.italic,
    #          label = paste0(fam.code,
    #                         " (",
    #                         Partition_,
    #                         ")"
    #                         )
    #         ) +
  scale_x_continuous(breaks = seq(0, max.clubsize, 20),
                     limits = c(0, max.clubsize)) +
  scale_y_continuous(breaks = seq(0, 4.5, .5),
                     limits = c(0, 4.2))

return(rc.plot)

}


font.cmu.serif.upright.italic <- "CMU Serif Upright Italic"
font.cmu.classical.serif <- "CMU Classical Serif"

col.title.height <- 0.2

plot.rc.whole <- plot_grid(
  plot_labels("  Whole"),
  plot_richclub(fam.code = "HC", Partition_ = "whole"),
  plot_richclub(fam.code = "HD", Partition_ = "whole"),
  plot_richclub(fam.code = "BL", Partition_ = "whole"),
  plot_richclub(fam.code = "SD", Partition_ = "whole"),
  plot_richclub(fam.code = "SC", Partition_ = "whole"),
  rel_heights = c(col.title.height,1,1,1,1,1),
  nrow = 6
)

plot.rc.majority <- plot_grid(
  plot_labels("  Majority"),
  plot_richclub(fam.code = "HC", Partition_ = "majority"),
  plot_richclub(fam.code = "HD", Partition_ = "majority"),
  plot_richclub(fam.code = "BL", Partition_ = "majority"),
  plot_richclub(fam.code = "SD", Partition_ = "majority"),
  plot_richclub(fam.code = "SC", Partition_ = "majority"),
  rel_heights = c(col.title.height,1,1,1,1,1),
  nrow = 6
)

plot.rc.minority <- plot_grid(
  plot_labels("  Minority"),
  plot_richclub(fam.code = "HC", Partition_ = "minority"),
  plot_richclub(fam.code = "HD", Partition_ = "minority"),
  plot_richclub(fam.code = "BL", Partition_ = "minority"),
  plot_richclub(fam.code = "SD", Partition_ = "minority"),
  plot_richclub(fam.code = "SC", Partition_ = "minority"),
  rel_heights = c(col.title.height,1,1,1,1,1),
  nrow = 6
)

plot.vertical.labels <- plot_grid(plot_labels("  ", angle = 90),
                                  plot_labels("HC", angle = 90),
                                  plot_labels("HD", angle = 90),
                                  plot_labels("BL", angle = 90),
                                  plot_labels("SD", angle = 90),
                                  plot_labels("SC", angle = 90),
                                  plot_labels("  ", angle = 90),
                                  rel_heights = c(0.5*col.title.height,
                                                  1,1,1,1,1,
                                                  0.25*col.title.height),
                                  nrow = 7
                                  )

plot.rc.all <- plot_grid(
          plot.vertical.labels,
          plot.rc.whole,
          plot.rc.majority,
          plot.rc.minority,
          ncol = 4,
          # hjust = -0.5,
          # vjust = 1.5,
          rel_widths = c(1,8.5,6.5,4))

plot.rc.final <- plot_grid(plot_labels("Normalized Rich-Club Coefficient",
                                    size = 7.5),
          plot.all,
          nrow = 2,
          rel_heights = c(1.5*col.title.height, 5)
          )

save_plot("Normalized-rich-club.pdf",
          plot.final,
          base_height = 11.69,
          base_width = 8.27)


tikz(file = "plot_test.tex", width = 8.27, height = 11.69)

print(plot.rc.all)

dev.off()


# Plotting connectivities -------------------------------------------------



# Trash pad ---------------------------------------------------------------


## Redoing family evolution plots, this time without RC/efficiency ---------

all.partitions <- snp$Partition %>% 
  unique() %>% 
  sort() %>% 
  as.character() %>%
  as.list()

all.vd <- snp$Verbal.Description %>% 
  unique() %>% 
  sort() %>% 
  as.character() %>%
  as.list()


system.time(
  for(vd in all.vd){
    for(part in all.partitions){
      extract_plotcoefs.glued(this.Verbal.Description = vd,
                              this.Partition = part, 
                              snp=snp,
                              coef.range = coefindex)
    }
  }
)


myFont <- "CMU Classical Serif"

ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point() +
  ggtitle("Fuel Efficiency of 32 Cars") +
  xlab("Weight (x1000 lb)") + ylab("Miles per Gallon") +
  theme(text = element_text(size = 16, family = myFont)) +
  annotate("text", 4, 30, label = 'Palatino Linotype',
           family = myFont, size = 10) +
  annotate("label", 1, 35,
           label = 'Roboto',
           hjust = 0,
           # vjust = 1,
           family = myFont, size = 10)

## On-screen device
print(a)

