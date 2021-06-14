
# Notes -------------------------------------------------------------------
#############################
######## 2021-06-14 #########
#############################
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


# Reading snp

load("./functions/functions_save_vars.R")

system.time(
  load("./data/snp_all_5k_20190712_1730.RData")
)

snp.tmp <- snp

# Correcting family labels and names

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

snp.tmp <- snp.tmp %>% 
  add_column(fn = paste0(rep(c(1:50),
                             each=800),
                         "_",
                         snp.tmp$Verbal.Description,
                         "_",
                         snp.tmp$Owner),
             .after = 6) %>%
  arrange(fn) %>%
  add_column(ID = paste0(snp.tmp$Verbal.Description %>%
                           as.character() %>% 
                           sort(decreasing = FALSE) %>%
                           substr(1,2),
                         rep(
                           rep(c(1:10),
                               each=800),
                           5)
                         ),
             .after = 1)
  
# Correcting small-world index (changing C*E to C/APL)

snp.tmp$Small.World <- snp.tmp$Clustering/snp.tmp$Average.Path.Length


# saving snp and variants back

snp <- snp.tmp
save_vars("snp", prefix = "snp_all_5k")

snp.lean <- snp.tmp %>% select(-adj.mat.vect)
save_vars("snp.lean", prefix = "snp.lean_all_5k")

snp <- snp.tmp %>% filter(Rewiring == 1000000)
save_vars("snp", prefix = "snp_only-1e+6")
