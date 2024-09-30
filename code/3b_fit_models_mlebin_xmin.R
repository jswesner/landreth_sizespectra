# Load the packages
library(tidyverse)
library(brms)
library(isdbayes)

# load data

landreth_fishmacros_data = readRDS(file = "data/landreth_fishmacros_data.rds") 

# refit models
brm_landreth_fishmacros_intercept_summin <- update(readRDS("models/brm_landreth_fishmacros_intercept.rds"),
                                                   newdata = landreth_fishmacros_data)
saveRDS(brm_landreth_fishmacros_intercept_summin, file = "models/brm_landreth_fishmacros_intercept_summin.rds")
brm_landreth_fishmacros_agtopo_summin <- update(readRDS("models/brm_landreth_fishmacros_agtopo.rds"),
                                                newdata = landreth_fishmacros_data)
saveRDS(brm_landreth_fishmacros_agtopo_summin, file = "models/brm_landreth_fishmacros_agtopo_summin.rds")
brm_landreth_fishmacros_trophic_summin <- update(readRDS("models/brm_landreth_fishmacros_trophic.rds"),
                                                 newdata = landreth_fishmacros_data)
saveRDS(brm_landreth_fishmacros_trophic_summin, file = "models/brm_landreth_fishmacros_trophic_summin.rds")

