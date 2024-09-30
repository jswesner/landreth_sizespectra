# Load the packages
library(tidyverse)
library(brms)
library(isdbayes)

# load data

landreth_fish_data = readRDS(file = "data/landreth_fish_data.rds") 

# refit models
brm_landreth_fish_intercept_summin <- update(readRDS("models/brm_landreth_fish_intercept.rds"),
                                                   newdata = landreth_fish_data, 
                                               cores = 4)
saveRDS(brm_landreth_fish_intercept_summin, file = "models/brm_landreth_fish_intercept_summin.rds")
brm_landreth_fish_agtopo_summin <- update(readRDS("models/brm_landreth_fish_agtopo.rds"),
                                                newdata = landreth_fish_data, 
                                            cores = 4)
saveRDS(brm_landreth_fish_agtopo_summin, file = "models/brm_landreth_fish_agtopo_summin.rds")
brm_landreth_fish_trophic_summin <- update(readRDS("models/brm_landreth_fish_trophic.rds"),
                                                 newdata = landreth_fish_data, 
                                             cores = 4)
saveRDS(brm_landreth_fish_trophic_summin, file = "models/brm_landreth_fish_trophic_summin.rds")

