# Klinkovská et al. (2024) Biological Conservation: Dynamics of the Czech flora over the last 60 years: winners, losers and causes of changes

## Supplementary code and data to the article: 

Klinkovská K., Glaser M., Danihelka J., Kaplan Z., Knollová I., Novotný P., Pyšek P., Řezníčková M., Wild J. & Chytrý M. (2024) Dynamics of the Czech flora over the last 60 years: winners, losers and causes of changes. Biological Conservation.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10625478.svg)](https://doi.org/10.5281/zenodo.10625478)

This repository contains data on plant species occurrences and plant characteristics from the Pladias Database of the Czech Flora and Vegetation (https://pladias.cz/en) and code used to analyse temporal trends of species of Czech flora from 1961 to 2020 using occupancy models. 

## Data
### Data used to calculate temporal trends of change in species frequency
The data used to estimate the trends of change in species frequency using occupancy models are available in the Zenodo repository ([![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10625478.svg)](https://doi.org/10.5281/zenodo.10625478)) and include a CSV file for each species, for which we calculated the occupancy model. The species name is in the file name. Each file consists of six comma-delimited columns: 
* `grid_year_auth` = visit ID, a combination of a grid cell, year, author and record origin category (converted to a numeric variable)
* `Nr_of_spec` = number of species recorded during the given visit 
* `half_dec` = half-decade (1 = 1960s, 2 = 1965s...)
* `grid_small` = grid cell ID (converted to numeric variable)
* `hab.map` = record origin category (1 = Natura 2000 mapping, 0 = other record)
* `value` = presence (1) or absence (0) of the given species

### Data used to compare differences in species characteristics for groups of species with different temporal trends
* `records_decades.csv`: Numbers of occurrences of each species in each half-decade.
* `occ_results.csv`: Results from the occupancy models summarized in one csv file for further analysis of differences between groups of species with different temporal trends.
* `not_conv_all.csv`: Appendix A.2: List of species for which the occupancy models did not converge.
* `species_characteristics.csv`: Species characteristics used to analyse the differences between groups of species with different temporal trends.
* `clustering_results.csv`: Assignment of species to a cluster according to the time series clustering algorithm.

## Code
### Occupancy model code for JAGS:
* `randomwalk_ht.txt`: Original code from Outhwaite et al. (2018, https://doi.org/10.1016/j.ecolind.2018.05.010) with the addition of record origin effect.

### R scripts:

* `occupancy_model.R`: Run occupancy model for the given species from R.
* `diag.R`: Extract diagnostics for occupancy models.
* `occ_results.R`: Get csv file with occupancy estimates from the RData object.
* `meta_results.R`: Subsequent analysis of species trends linked with species characteristics - linear models, chi-square tests, time series clustering.
