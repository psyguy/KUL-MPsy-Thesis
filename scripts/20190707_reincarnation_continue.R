#!/usr/bin/env Rscript

Args <- commandArgs(TRUE)
index <- as.numeric(Args[1])

# index <- 7
# first making a brain, from 0626_hpc.R -----------------------------------

source("./functions/functions_reports.R")

## read the HPC instructions here:
## https://github.com/psyguy/Emotion-Dynamics-1/blob/master/ED%201%20-%20Codes/3.%20mirt-Model-Comp/correct_HPC-ready_May22/HPC%20readme.txt

sampled.path <- "./data/5200-edges/"

## used locally to manually make the following vector
# sampled.names <-
#   list.files(path = sampled.path, pattern = "*life-01")
# 
# names <- sampled.names %>%
#   my_gsub(".*-5.2k_") %>%
#   my_gsub("_.*")

names <- c("Sam Evrard", "Jordan Vermeersch",
           "Dorian Gerard", "Daan Pauwels",
           "Arthur Lejeune", "Sam Verheyen",
           "Alexandre Laurent", "Maxim Declercq",
           "Loic Van De Velde", "Lauren Mathieu",
           "Ines Verheyen", "Jana Wouters",
           "Celine Thys", "Thibault Janssens",
           "Lea Vandenberghe", "Jordan De Pauw",
           "Emilie Petit", "Samuel Devos",
           "Anna Dupont", "Elise De Cock",
           "William De Meyer", "Maxime Simon",
           "Romain Thys", "Alexandre Van Dyck",
           "Oceane Bauwens", "Yana Van Damme",
           "Kevin Jansen", "Milan Jacobs",
           "Morgane Leclercq", "Julien Thomas",
           "Noah Wuyts", "Logan Wauters",
           "Milan Verlinden", "Loic Vermeiren",
           "Oceane Pauwels", "Ruben Declercq",
           "Niels De Cock", "Anouk De Pauw",
           "David Cornelis", "Martin Devos",
           "Ines Lemaire", "Celine Simon",
           "Seppe Jansen", "Ine Verheyen",
           "Margaux Bogaert", "Hannah Bauwens")


name <- names[index]
t0 <- Sys.time()

for (reincarnation in 8:10) {
  t1 <- Sys.time()
  
  pattern <- "_g-0.3k-5.2k"
  
  sampled.names <-
    list.files(path = sampled.path, pattern = "*.RData")
  
  # names <- sampled.names %>%
  #   my_gsub(paste0(".*",
  #                  pat.tmp,
  #                  "_"), "") %>%
  #   my_gsub("\\_.*", "") %>%
  #   sort()
  
  this.owner.oldest <-
    sampled.names[grepl(name, sampled.names)] %>%
    sort() %>%
    tail(1)
  
  current.life <-
    this.owner.oldest %>% substring(6, 7) %>% as.numeric()
  
  # for(sampled in r.this){
  
  life.prefix <- sprintf("life-%02d", (current.life + 1))
  
  load(paste(sampled.path, this.owner.oldest, sep = ""))
  
  ## rebirth
  brain_case@initial <- brain_case@now
  # making initial$mat.con a list (standard form)
  brain_case@initial$mat.connectivity <-
    brain_case@now$mat.connectivity %>% list()
  # removing the long history of previous life
  brain_case@history <- brain_case@initial
  
  brain_case <- partition_culture(brain_case = brain_case,
                                  final.age = (brain_case@age$rewires %/% 1000) +
                                    100)
  
  Sys.time() - t1
  
  # tryCatch({
  #   reports_netviz(brain_case)
  # }, error = function(e) {
  #   print(paste("Error plotting", this.owner.oldest))
  # })
  #
  # Sys.time() - t1
  
  
  save_vars(
    list.of.vars = "brain_case",
    prefix = paste(
      life.prefix,
      brain_case@parameters$brain.code,
      brain_case@name,
      sep = "_"
    ),
    path = sampled.path
  )
  
}

Sys.time() - t0
    