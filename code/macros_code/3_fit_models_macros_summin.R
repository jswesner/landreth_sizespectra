# Load the packages
library(tidyverse)
library(brms)
library(isdbayes)

# load data

peak_min_sizes <- readRDS("data/peak_min_sizes_macros.RDS")
landreth_macros_data = readRDS(file = "data/landreth_macros_data.rds") %>% 
  left_join(peak_min_sizes) %>% 
  filter(dw_g >= sum_min) %>% 
  group_by(stream) %>% 
  mutate(xmin = min(dw_g),
         xmax = max(dw_g))

saveRDS(landreth_macros_data, file = "data/landreth_macros_data.rds")

# refit models
brm_landreth_macros_intercept_summin <- update(readRDS("models/brm_landreth_macros_intercept.rds"),
                                                   newdata = landreth_macros_data, 
                                               cores = 4)
saveRDS(brm_landreth_macros_intercept_summin, file = "models/brm_landreth_macros_intercept_summin.rds")
brm_landreth_macros_agtopo_summin <- update(readRDS("models/brm_landreth_macros_agtopo.rds"),
                                                newdata = landreth_macros_data, 
                                            cores = 4)
saveRDS(brm_landreth_macros_agtopo_summin, file = "models/brm_landreth_macros_agtopo_summin.rds")
brm_landreth_macros_trophic_summin <- update(readRDS("models/brm_landreth_macros_trophic.rds"),
                                                 newdata = landreth_macros_data, 
                                             cores = 4)
saveRDS(brm_landreth_macros_trophic_summin, file = "models/brm_landreth_macros_trophic_summin.rds")

