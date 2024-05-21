library(tidyverse)
library(brms)
library(isdbayes)

library(tidyverse)
library(brms)
library(isdbayes)

# model fish only ---------------------------------------------------------
# load data
landreth_fishmacros_data = read_csv(file = "data/landreth_fishmacros_data.csv") 

brm_landreth_fishmacros_topo = brm(dw_g|vreal(counts, xmin, xmax) ~ p_elva_std_s + p_mda_slope_s + p_da_s + (1|sample_id),
                              data = landreth_fishmacros_data,
                              stanvars = stanvars,
                              family = paretocounts(),
                              prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
                                        prior(exponential(5), class = "sd"),
                                        prior(normal(0, 1), class = "b")),
                              iter = 2000, chains = 4)

saveRDS(brm_landreth_fishmacros_topo, file = "models/brm_landreth_fishmacros_topo.rds")