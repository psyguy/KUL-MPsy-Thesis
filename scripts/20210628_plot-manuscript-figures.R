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

