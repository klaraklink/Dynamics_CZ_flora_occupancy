#' Supplementary code to the article: 
#' Klinkovska et al. (2024) Dynamics of the Czech flora over the last 60 years: winners, losers and causes of changes. Biological Conservation.
#' 
#' Author: Michael Glaser, Klara Klinkovska, 2023-06-21
#' R version 4.0.0
#' run at MetaCentrum 

library(readr)
library(backports)
library(lattice)
library(runjags)

args = commandArgs(trailingOnly = TRUE)

### function for occupancy model
occ.model <- function(n.iter, n.chains, n.burnin, n.thin, n.adapt, save.pars, specname){

  sysdate <- gsub("-","",Sys.Date())

  # read csv file for one species
  specdata <- read_csv(paste0("species/", specname, ".csv"))
  
  ### create vectors for data.list
  nyear       <- length(unique(specdata$half_dec))            ### number of years 
  nsite       <- max(specdata$grid_small)                       ### number of sites 
  nvisit      <- max(specdata$grid_year_auth)                ### number of visits 
  y           <- specdata$value                       ### detection status of visit
  logL        <- log(specdata$Nr_of_spec)                   ### logarithm of list length 
  Site        <- specdata$grid_small                            ### Site associated with visit
  Year        <- specdata$half_dec                            ### Year associated with visit
  h           <- specdata$hab.map                           ### record type, habitat mapping or not
  
  ### combine into list (JAGS needs a list as input format)
  data.list <- list("nyear" = nyear, "nsite" = nsite, "nvisit" = nvisit,
                  "y" = y, "logL" = logL, "Site" = Site, "Year" = Year, "h" = h)
  
  ### Occupancy model
  mod.start <- Sys.time()
  
  jagsres <- autorun.jags(model = paste0("MODELS/randomwalk_ht.txt"), method = "parallel",
                          max.time = "94 hours", data = data.list, modules = c("glm","dic"), 
                          monitor = save.pars, thin = n.thin, thin.sample = T, adapt = n.adapt, 
                          startburnin = n.burnin, n.chains = n.chains, startsample = n.iter)
  
  mod.end <- Sys.time()
  
  runtime <- as.numeric(mod.end-mod.start, units = "hours")
  
  saveRDS(jagsres, file = paste0("OUTPUT/", specname, "_pladias_randomwalk_grid_small_thin", n.thin,".rds"))
  
  write_csv(as.data.frame(runtime), paste0("OUTPUT/", specname, "_runtime.csv"))
  
}

# model parameters
n.iter    <- 4000  
n.chains  <- 2 
n.thin    <- 5 
n.adapt   <- 500 
n.burnin  <- 5000 
save.pars <- "psi.fs"
specname <- args[1]

# run the function
thisspec.occ <- occ.model(n.iter, n.chains, n.burnin, n.thin, n.adapt, save.pars, specname)