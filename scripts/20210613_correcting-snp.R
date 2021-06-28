
# Notes -------------------------------------------------------------------
##########################################################
######## 2021-06-14 ######################################
##########################################################
##
## snp data objects prior to 2021 had two flaws:
## 1.family names must be corrected:
##    - Under coupled UC --> sub-coupled SC
##    - Over coupled OC --> hyper-coupled HC
##    - over-turbulent OT --> hyper-divergent HD
##    - underturbulent UT --> sub-divergent SD
## 2. small-world index is calculated incorrectly (changing C*E to C/APL)
##
## So we first load the complete snp, make changes, and save versions of it
## These will take time (~10 minutes)
##
##########################################################
######## 2021-06-14 ######################################
##########################################################
##
## The ordering of brains (1-50) and families is corrected to match 2019 files.
## So use the following from now on:
##    - snp_all_5k_20210628_1238.RData (zipped into three parts, for GitHub)
##    - snp.lean_all_5k_20210628_1241.RData
##    - snp_only-1e+6_20210628_1241.RData
## 
##########################################################

# Reading snp -------------------------------------------------------------

source("./functions/functions_save_vars.R")

system.time(
  load("./data/snp_all_5k_20190712_1730.RData")
)

snp.tmp <- snp


# Correcting family labels and names --------------------------------------

snp.tmp$Verbal.Description <-
  plyr::mapvalues(snp.tmp$Verbal.Description,
                 from = c("Hyper-coupled minority",
                           "Hyper-chaotic minority",
                           "Homogeneous Society",
                           "Hypo-chaotic minority",
                           "Hypo-coupled minority"),
                  to = c("HC_Hyper-coupled minority",
                         "HD_Hyper-divergent minority",
                         "BL_Baseline homogeneous",
                         "SD_Sub-divergent minority",
                         "SC_Sub-coupled minority")
                 )

# To sort everything alphabetically (instead of HC->HD->BL->SD->SC sorting)

snp.tmp$Verbal.Description <- snp.tmp$Verbal.Description %>%
  as.character()
snp.tmp$Owner <- snp.tmp$Owner %>% as.character()
# Note that the abovel lines changes plot orders in the exported PDFs,
# which is corrected a few lines below.


# Making family codes and individual IDs ----------------------------------

snp.tmp <- snp.tmp %>% 
  arrange(Verbal.Description, Owner) %>% 
  mutate(fn = paste0(
                     str_pad(# this is to "correct" ordering for PDFs
                             rep(c(21:30,
                                   1:10,
                                   11:20,
                                   41:50,
                                   31:40),
                                 each = 800),
                             2,
                             pad = "0"),
                     "_",
                     Verbal.Description,
                     "_",
                     Owner),
             .after = 6) %>%
  arrange(fn) %>%
  mutate(ID = paste0(Verbal.Description %>%
                     as.character() %>% 
                     substr(1,2),
                     rep(rep(c(1:10),
                             each = 800),
                         5)
                     ),
         .after = 1)


# Correcting small-world index (changing C*E to C/APL) --------------------

snp.tmp$Small.World <- snp.tmp$Clustering/snp.tmp$Average.Path.Length


# saving snp and variants back --------------------------------------------

snp <- snp.tmp
save_vars("snp", prefix = "snp_all_5k")

snp.lean <- snp.tmp %>% select(-adj.mat.vect)
save_vars("snp.lean", prefix = "snp.lean_all_5k")

snp <- snp.tmp %>% filter(Rewiring == 1000000)
save_vars("snp", prefix = "snp_only-1e+6")
