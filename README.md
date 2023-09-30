# Occupancy_models_CZ

Supplementary code to the article: 

Klinkovská et al. Dynamics of the Czech flora over the last 60 years: winners, losers and causes of changes

Author: Klára Klinkovská, Michael Glaser 

## R scripts:

1. data_preparation.R - loading data exported from Pladias, preparation for Occupancy models
2. occupancy_model.R - run occupancy model at the server (MetaCentrum)
3. diag.R - occupancy model diagnostics
4. occ_results.R - get csv with occupancy estimates from RData object
5. main.R - analysis of species trends related with their characteristics, linear models, chi-square tests, time series clustering

## Python scripts:
* jobcreate.py - create bash job for MetaCentrum

## Bash scripts:
1. diag_meta.sh - get diagnostics of occupancy models
2. rdata_to_csv - get csv with occupancy estimates from RData object


