get_networks <- function(networks,
                         organism = c("mouse", "human"),
                         top_edges_only = TRUE) {
  
  if(!dir.exists(paste(here("data/"), paste(organism, "networks", sep = "_"), sep = "/"))) {
    
    dir.create(paste(here("data/"), paste(organism, "networks", sep = "_"), sep = "/"))
    
    cat(paste(organism, "network directory created, checking for networks to be dowloaded...\n", sep = " "))
    
  }
  
  cat(paste(organism, "-specific directory exists, checking for networks to be downloaded...\n", sep = ""))
  
  nets_to_download <- vector() 
  
  for(i in 1:length(networks)) {
    
    if(!file.exists(paste(networks[i], "top.gz", sep = "_"))) {
      
      nets_to_download[i] <- networks[i]
      
    }
    
  }
  
  cat(paste(length(nets_to_download), organism, "networks will be downloaded...\n"))
  
  for(i in 1:length(nets_to_download)) {
    
    download.tissue.net(tissue = nets_to_download[i], 
                        organism = organism,
                        top.edges.only = top_edges_only,
                        project.dir = paste(here("data/"), paste(organism, "networks", sep = "_"), sep = "/"))
    
  }
  
  cat("Networks have been downloaded successfully!")
  
}

