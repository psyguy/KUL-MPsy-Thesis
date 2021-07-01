# Notes -------------------------------------------------------------------
##########################################################
######## plotting adj+net with PCA_angle seriation #######
######## 2021-06-14 ######################################
##########################################################
## 
## From commit 2aa903a on, the following snp data files should be used:
##    - snp_all_5k_20210628_1238.RData (zipped into three parts, for GitHub)
##    - snp.lean_all_5k_20210628_1241.RData
##    - snp_only-1e+6_20210628_1241.RData
## 
##########################################################


# Loading functions and snp -----------------------------------------------

library(librarian)
shelf(cowplot)
source("functions/functions_my.R")
source("functions/functions_netmeas.R")
source("./functions/functions_extract.R")
load("./data/snp_only-1e+6_20210628_1241.RData")
load("./data/snp.lean_all_5k_20210628_1241.RData")

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


save_vars("snp.rc", prefix = "Lean snp at 1e+6 with rich club")


# Plotting rich club ------------------------------------------------------


make.df <- function(input.df, col){
  
  d <- input.df %>% 
    select("Partition",
           "ID",
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
                           "ID",
                           "Verbal.Description",
                           col,
                           "Club Size")
  output.df %>% return()
}


# plotting RC values and significance levels ------------------------------

rc.res <- snp.rc

# rc.res$num <- rc.res$num %>% as.character() %>% as.factor()

plot_richclub <- function(rc.res = snp.rc,
                          fam.code = "BL",
                          Partition_ = "whole"){
  
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


rc.plot <- rc.this %>%  ggplot(
  aes(x = `Club Size`,
      y = `Normalized Rich Club`,
      colour = ID)) +
  geom_line(size = 2, alpha = 0.35) +
  # ggplot2::xlim(0,75) +
  # ggplot2::ylim(0.5,2.20) +
  theme_bw() +
  ggtitle(paste("Normalized Rich Club for", fam.code)) +
  theme(plot.title = element_text(hjust = 0.5),
        # axis.title.y = element_blank(),
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        # axis.line = element_line(colour = "black"),
        # panel.background = element_blank(),
        legend.position = "none"
        ) +
  geom_point(aes(x = `Club Size`,
                 y = rc.p,
                 color = ID),
             data = rc.this,
             size = 1.2,
             shape = 16) +
  geom_hline(yintercept = 1,linetype="dashed") +
  scale_x_continuous(breaks = seq(0, max.clubsize, 20),
                     limits = c(0, max.clubsize)) +
  scale_y_continuous(breaks = seq(0, 4.5, .5),
                     limits = c(0, 4.2))

return(rc.plot)

}

plot_richclub(fam.code = "BL")
# p.bl <- 
plot_richclub(fam.code = "HD", Partition_ = "majority")

plot_grid(p.sd,p.bl)



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

a <- ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point() +
  ggtitle("Fuel Efficiency of 32 Cars") +
  xlab("Weight (x1000 lb)") + ylab("Miles per Gallon") +
  theme(text = element_text(size = 16, family = myFont)) +
  annotate("text", 4, 30, label = 'Palatino Linotype',
           family = myFont, size = 10) +
  annotate("text", 1, 11, label = 'Roboto', hjust = 0,
           family = myFont, size = 10)

## On-screen device
print(a)

