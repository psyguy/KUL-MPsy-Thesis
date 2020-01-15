# remotes::install_github("centerforopenscience/osfr")
library(osfr)
osf_auth("L6b9zDbSgmQGGVyQMyXc73hr1lWHQKWtlgTunYVSGgDEMrCqDfn3pAnx3zv0qc9xGpUHXM")
osf_retrieve_node("625d8") %>%
  osf_ls_files(path = "Figures/Network maturation",
               n_max = 50) %>%
  osf_download()

l.f <- 

pdf.names <- list.files(pattern = "Network maturation")
d <- pdf.names %>% 
  gsub("Network maturation of ", "", .) %>%
  gsub(" \\(.*","",.) %>%
  as.data.frame()
d <- pdf.names %>%
  gsub(".* \\(","",.) %>% 
  gsub("\\).pdf","",.) %>%
  cbind(d, pdf.names)

colnames(d) <- c("family", "name", "pdf.names")

dd <- d %>%
  mutate(fn = paste(family,name,sep="_")) %>%
  arrange(fn) %>%
  mutate(fn = paste0(rep(c(3,2,1,4,5),each=10), "_",fn,".pdf"))

d_ply(dd[1:2,],1,function(x){select_pages(40, x$pdf.names, x$fn)})

for(i in 1:50){
  select_pages(40, as.character(dd[i,3]), as.character(dd[i,4]))
}

staple_pdf(input_files = sort(dd$fn), output_filepath = "networks at 1e+6.pdf")

file.rename(list.files(path = "./figures/rc"), paste0("./figures/rc/",dd,".pdf"))

file.rename(list.files(pattern="water_*.img"), paste0("water_", 1:700))


d %>%
  mutate(filename = paste0("Final Rich Club and Betweenness of ",
                           name,
                           " (",
                           family)
  ) %>% 
  pull(filename) #%>% 
paste0("./figures/",.) %>% 
  staplr::staple_pdf(output_filepath = paste("meh")
  )

