# Klinkovská et al. (2024) Biological Conservation: Dynamics of the Czech flora over the last 60 years: winners, losers and causes of changes

## Supplementary code to the article: 

Klinkovská K., Glaser M., Danihelka J., Kaplan Z., Knollová I., Novotný P., Pyšek P., Řezníčková M., Wild J. & Chytrý M. (2024) Dynamics of the Czech flora over the last 60 years: winners, losers and causes of changes. Biological Conservation.

Author: Klára Klinkovská (klinkovska.klara@gmail.com), Michael Glaser (michael.glaser@univie.ac.at)

This repository contains R scripts used to analyse temporal trends in the data from the Pladias Database of the Czech Flora and Vegetation (https://pladias.cz/en/).

## R scripts:

* `occupancy_model.R`: Calculation of occupancy model for the given species.
* `diag.R`: Extract diagnostics for occupancy models.
* `occ_results.R` Get csv file with occupancy estimates from the RData object.
* `meta_results.R` Subsequent analysis of species trends linked with species characteristics - linear models, chi-square tests, time series clustering.

model for JAGS - randomwalk txt file
