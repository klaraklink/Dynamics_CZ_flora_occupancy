#' Supplementary code to the article: 
#' Klinkovska et al. (2024) Dynamics of the Czech flora over the last 60 years: winners, losers and causes of changes. Biological Conservation.
#' 
#' Author: Michael Glaser, Klara Klinkovska, 2023-06-21
#' R version 4.0.0
#' run at MetaCentrum

library(tidyverse)
library(ggmcmc)
library(backports)
library(grDevices)
library(ggforce)

# function to get diagnostics
occ.diag <- function(filepath){
  
  RData.list <- list.files(path = filepath,pattern = paste0("_pladias_randomwalk_grid_small_thin5.rds"))
  
  conv.df <- data.frame(Parameter = character(),
                        Chain = numeric(),
                        Diagnostic = character(),
                        value = numeric(),
                        specname = character())
  
  
  jagsres.ggs <- data.frame(Iteration = numeric(), Chain = numeric(),
                            Parameter = character(), value = numeric(),
                            ParameterOriginal = character())
  
  for (e in 1:length(RData.list)){
    
    ### load RData object
    jagsres <- readRDS(paste0(filepath, "/", RData.list[e]))
    specname <- str_remove(RData.list[e], "_pladias_randomwalk_grid_small_thin5.rds") %>% 
      str_replace_all("_", " ")
    parlabeldf <- data.frame(Parameter = paste0("psi.fs[", 1:12, "]"),
                             Label = paste0(seq(1960, 2015, 5), "s"))
    
    ### extract diagnostics
    jagsres.ggs2 <- ggs(jagsres$mcmc,keep_original_order = T, family = "psi.fs", par_labels = parlabeldf)
    diagn <- ggs_diagnostics(jagsres.ggs2)
    diagn$specname <- specname
    conv.df <- rbind(diagn, conv.df)
    jagsres.ggs2$specname <- specname
    jagsres.ggs <- rbind(jagsres.ggs, jagsres.ggs2)
  }
  
  tbl.list <- list(conv = conv.df, jagsres.ggs = jagsres.ggs
  )
  
  return(tbl.list)
  
}

# get diagnostics
test <- occ.diag("OUTPUT/")

# list of species with not converging results according to Rhat statistics
n_conv <- test$conv %>% 
  filter(Diagnostic == "Rhat" & value >= 1.1) %>% 
  write_csv("diag/not_conv2.csv")

# export plots to check diagnostics manually
specunique <- unique(test$conv$specname)
for (sp in specunique){
  
  jags.ggs.sp <- test$jagsres.ggs[test$jagsres.ggs$specname == sp,]
  
  pdf(paste0("diag/", str_replace_all(sp, " ", "_"), "_pladias_randomwalk_grid_small_thin5.pdf")) 
  
  # traceplots
  for (pages in c(1:5)){
    thispage <- ggplot(jags.ggs.sp, aes(x = Iteration, y = value, color = factor(Chain)))+
      theme_classic()+
      geom_line() +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1))+
      scale_x_continuous(expand = c(0, 0))+
      facet_wrap_paginate(~Parameter, nrow = 2, ncol = 1, scales = "fixed", page = pages)+
      ggtitle(sp)
    
    print(thispage)
  }
  
  print(ggs_geweke(jags.ggs.sp)+
          ggtitle(paste("Geweke Plot for", sp)))
  
  # autocorrelation plots
  print(ggs_autocorrelation(jags.ggs.sp))+
    ggtitle(sp)
  
  # occupancy estimates, CIs and linear trend
  jagsres.ci <- ci(test$jagsres.ggs[test$jagsres.ggs$specname == sp,])
  jagsres.ci$parnum <- as.numeric(as.character(gsub("s", "", jagsres.ci$Parameter)))
  
  occ.lm <- lm(jagsres.ci$median~jagsres.ci$parnum)
  occ.int <- occ.lm$coefficients[["(Intercept)"]]
  occ.slo <- occ.lm$coefficients[["jagsres.ci$parnum"]]
  
  occ.res <- ggplot(jagsres.ci, aes(x = parnum, y = median,))+ theme_classic()+
    geom_point(size = 2.5, shape = 15)+
    geom_abline(slope = occ.slo, intercept = occ.int, col = "red", size = 2)+
    geom_line(linetype = "dashed")+
    geom_ribbon(aes(ymin = Low,ymax = High), color = "grey20", alpha = 0.1)+
    geom_ribbon(aes(ymin = low,ymax = high), color = "grey50", alpha = 0.1)+
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.1))+
    scale_x_continuous(limits = c(min(jagsres.ci$parnum),max(jagsres.ci$parnum)),expand = c(0, 0.5),
                       breaks = jagsres.ci$parnum, labels = jagsres.ci$Parameter)+
    theme(axis.text.x = element_text(hjust = 0.85))+
    ggtitle(paste("Occupancy and CIs for", sp, "lm-trend = ", round(occ.slo, 3)))+
    xlab("timestep")+
    ylab("median occupancy")
  
  print(occ.res)
  
  dev.off()
  
}
