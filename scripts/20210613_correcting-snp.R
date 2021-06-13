
# Notes -------------------------------------------------------------------

## snp data has two flaws:
## 1. small-world index is calculated incorrectly (changing C*E to C/APL)
## 2.family names must be corrected:
##    - Under coupled UC --> sub-coupled SC
##    - Over coupled OC --> hyper-coupled HC
##    - over-turbulent OT --> hyper-divergent HD
##    - underturbulent UT --> sub-divergent SD

## So we first load the complete snp, make changes, and save versions of it
# system.time(
#   load("./data/snp_all_5k_20190712_1730.RData")
# )

snp.tmp <- snp

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
                  
(family.code <- snp.tmp$Verbal.Description %>%
  as.character() %>% 
  sort(decreasing = FALSE) %>% 
  substr(1,2) %>% 
  rep(each=10) %>% 
  paste0(rep(rep(c(1:50), each=800)),
         "_",
         .)
  )

(uu <- family.code %>% unique())

rep()

snp.tmp2 <- snp.tmp %>% 
  add_column(fn = paste0(snp.tmp$Verbal.Description,
                         "_",
                         snp.tmp$Owner),
             .after = 6) %>%
  arrange(fn) %>%
  # select(-adj.mat.vect) %>% 
  add_column(fn = paste0(rep(c(3,2,1,4,5),each=10),
                         "_",
                         snp.tmp$fn),
             .after = 7)
  arrange(fn) %>%
  cbind(family.codes = families, num =1:50) %>%
  mutate(newcodes = paste0(num,"_",family.codes,rep(c(1:10), 5), "_", Owner)) %>%
  mutate(famcode.num = paste0(family.codes, rep(c(1:10), 5)))


  
xx <- 1
ddd[(800*xx-8):(800*xx+8),]
  


# families <<- c("OC", "OT", "BL", "UT", "UC") %>% rep(each = 10)
families <<- c("HC", "HD", "BL", "SD", "SC") %>% rep(each = 10)

### correcting collumn names in snp

# correcting small-world index (changing C*E to C/APL)

snp.lite <- snp[sample(50,10),1:16]

snp$Small.World <- snp$Clustering/snp$Average.Path.Length
