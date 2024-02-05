#' Supplementary code to the article: 
#' Klinkovska et al. (2024) Dynamics of the Czech flora over the last 60 years: winners, losers and causes of changes. Biological Conservation.
#' 
#' Author: Michael Glaser, Klara Klinkovska, 2023-09-22
#' R version 4.0.0
#' ran at MetaCentrum

library(stringr)
library(ggmcmc)
library(readr)
library(dplyr)
library(data.table)

occ.res.occ <- function(filepath){
  
  RData.list <- list.files(path = filepath, pattern = ".rds")
  
  occ.df <- data.frame(Parameter = character(),
                     low = numeric(),
                     Low = numeric(),
                     median = numeric(),
                     High = numeric(),
                     high = numeric(),
                     ParameterOriginal = character(),
                     parnum = numeric(),
                     specname = character())
  
  for (e in 1:length(RData.list)){
    
    ### load RData object & extract specname
    jagsres <- readRDS(paste0(filepath, RData.list[e]))
    specname <- RData.list[e] %>% 
      str_remove("_pladias_randomwalk_grid_small_thin5.rds") %>% 
      str_replace_all("_", " ")
    
    parlabeldf <- data.frame(Parameter = paste0("psi.fs[", 1:12, "]"),
                           Label = paste0(seq(1960, 2015, 5), "s"))
    
    
    jagsres.ggs <- ggs(jagsres$mcmc, keep_original_order = T, family = "psi.fs", par_labels = parlabeldf)
    
    # calculate confidence intervals
    jagsres.ci <- ci(jagsres.ggs) 
    jagsres.ci$parnum <- as.numeric(as.character(gsub("s", "", jagsres.ci$Parameter)))
    
    ### extract occupancy estimates
    thisocc <- jagsres.ci
    thisocc$specname <- specname
    
    occ.df <- rbind(occ.df, thisocc)
    
    ### remove the RData object to unclutter environment
    rm(jagsres)
    
  }
  
  write_csv(occ.df, "results/occ_results.csv")
  
}

res.occ <- occ.res.occ("OUTPUT/")