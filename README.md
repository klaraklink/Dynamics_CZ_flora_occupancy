# Klinkovská et al. (2024) Biological Conservation: Dynamics of the Czech flora over the last 60 years: winners, losers and causes of changes

## Supplementary code and data to the article: 

Klinkovská K., Glaser M., Danihelka J., Kaplan Z., Knollová I., Novotný P., Pyšek P., Řezníčková M., Wild J. & Chytrý M. (2024) Dynamics of the Czech flora over the last 60 years: winners, losers and causes of changes. Biological Conservation.

This repository contains data from the Pladias Database of the Czech Flora and Vegetation (https://pladias.cz/en/) and code used to analyse temporal trends of species of Czech flora. 

## Data
The data include a CSV file with columns delimited by commas for each species, for which we calculated the occupancy model. The species name is included in the file name. Each file consists of six columns: 
* `grid_year_auth` = Visit ID, a combination of a grid cell, year, author and record origin category (converted to a numeric variable)
* `Nr_of_spec` = number of species recorded during the given visit 
* `half_dec` = half-decade (1 = 1960s, 2 = 1965s...)
* `grid_small` = grid cell ID (converted to numeric variable)
* `hab.map` = record origin category (1 = Natura 2000 mapping, 0 = other record)
* `value` = presence (1) or absence (0) of the given species

## Code
### Occupancy model code for JAGS:
* `randomwalk_ht.txt`: Original code from Outhwaite et al. (2018, https://doi.org/10.1016/j.ecolind.2018.05.010) with the addition of record origin effect

### R scripts:

* `occupancy_model.R`: Run occupancy model for the given species from R.
* `diag.R`: Extract diagnostics for occupancy models.
* `occ_results.R` Get csv file with occupancy estimates from the RData object.
* `meta_results.R` Subsequent analysis of species trends linked with species characteristics - linear models, chi-square tests, time series clustering.
