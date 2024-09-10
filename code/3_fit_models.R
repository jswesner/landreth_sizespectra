options(repos = c(CRAN = "https://cloud.r-project.org"))

# Check if 'tidyverse' is installed
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

# Check if 'brms' is installed
if (!requireNamespace("brms", quietly = TRUE)) {
  install.packages("brms")
}

# Check if 'remotes' is installed (needed to install from GitHub)
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Check if 'isdbayes' is installed, if not, install from GitHub
if (!requireNamespace("isdbayes", quietly = TRUE)) {
  remotes::install_github("jswesner/isdbayes")
}

# Load the packages
library(tidyverse)
library(brms)
library(isdbayes)

# load data
landreth_fish_data = readRDS(file = "data/landreth_fish_data.rds") %>% 
  mutate(counts = as.integer(counts*sample_area_m2))
# check for NAs in trophic groups
# landreth_fish_data %>% group_by(stream) %>% filter(is.na(fishmacros_invertivore_s))

# brm_landreth_macros_topo = brm(dw_g|vreal(counts, xmin, xmax) ~ p_mda_elev_s + p_mda_slope_s + p_da_s +
#                                  (1|stream) + (1|watershed),
#                                data = landreth_macros_data,
#                                stanvars = stanvars,
#                                family = paretocounts(),
#                                prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
#                                          prior(exponential(5), class = "sd"),
#                                          prior(normal(0, 1), class = "b")),
#                                iter = 2000, chains = 4)

# saveRDS(brm_landreth_macros_topo, file = "models/brm_landreth_macros_topo.rds")

brm_dummy = readRDS("models/brm_landreth_macros_topo.rds")

# brm_elev_std_fish = update(brm_dummy, formula = . ~ p_elva_std_s + (1|stream) + (1|watershed),
#                            newdata = landreth_fish_data,
#                            iter = 2000, chains = 4,
#                            prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
#                                      prior(exponential(5), class = "sd")),
#                            cores = 4)
# saveRDS(brm_elev_std_fish, file = "models/brm_elev_std_fish.rds")
# 
# brm_slope_std_fish = update(brm_dummy, formula = . ~ p_slope_std_s + p_mda_slope_s + (1|stream) + (1|watershed),
#                             newdata = landreth_fish_data,
#                             iter = 2000, chains = 4,
#                             prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
#                                       prior(exponential(5), class = "sd")),
#                             cores = 4)
# saveRDS(brm_slope_std_fish, file = "models/brm_slope_std_fish.rds")
# 
# brm_herbivore = update(brm_dummy, formula = . ~ fishmacros_herbivore_s + (1|stream) + (1|watershed),
#                        newdata = landreth_fish_data,
#                        iter = 2000, chains = 4,
#                        prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
#                                  prior(exponential(5), class = "sd")), 
#                        cores = 4)
# saveRDS(brm_herbivore, file = "models/brm_fish_herbivore.rds")
# 
# brm_invertivore = update(brm_dummy, formula = . ~ fishmacros_invertivore_s + (1|stream) + (1|watershed),
#                          newdata = landreth_fish_data,
#                          iter = 2000, chains = 4,
#                          prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
#                                    prior(exponential(5), class = "sd")), 
#                          cores = 4)
# saveRDS(brm_invertivore, file = "models/brm_fish_invertivore.rds")
# 
# brm_omnivore = update(brm_dummy, formula = . ~ fishmacros_omnivore_s + (1|stream) + (1|watershed),
#                       newdata = landreth_fish_data,
#                       iter = 2000, chains = 4,
#                       prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
#                                 prior(exponential(5), class = "sd")), 
#                       cores = 4)
# saveRDS(brm_omnivore, file = "models/brm_fish_omnivore.rds")
# 
# 
# brm_landreth_fish_intercept = update(brm_dummy,
#                                      newdata = landreth_fish_data, formula = . ~ 1 + (1|stream) + (1|watershed),
#                                      prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
#                                                prior(exponential(5), class = "sd")))
# 
# saveRDS(brm_landreth_fish_intercept, file = "models/brm_landreth_fish_intercept.rds")

# brm_landreth_fish_watershedpca1 = update(brm_dummy,
#                                          newdata = landreth_fish_data,
#                                          formula = .~ WatershedPC1_s + (1 | stream) + (1 | watershed),
#                                          prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
#                                                    prior(exponential(5), class = "sd")),
#                                          chains = 4, iter = 2000)
# 
# saveRDS(brm_landreth_fish_watershedpca1, file = "models/brm_landreth_fish_watershedpca1.rds")
# 
# brm_landreth_fish_watershedpca2 = update(brm_dummy,
#                                          newdata = landreth_fish_data,
#                                          formula = .~ WatershedPC2_s + (1 | stream) + (1 | watershed),
#                                          prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
#                                                    prior(exponential(5), class = "sd")),
#                                          chains = 4, iter = 2000)
# 
# saveRDS(brm_landreth_fish_watershedpca2, file = "models/brm_landreth_fish_watershedpca2.rds")
# 
# brm_landreth_fish_watershedpca12345 = update(brm_dummy,
#                                              newdata = landreth_fish_data,
#                                              formula = .~ WatershedPC1_s + WatershedPC2_s + WatershedPC3_s +
#                                                WatershedPC4_s + WatershedPC5_s + (1 | stream) + (1 | watershed),
#                                              prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
#                                                        prior(exponential(5), class = "sd")),
#                                              chains = 4, iter = 2000)
# 
# saveRDS(brm_landreth_fish_watershedpca12345, file = "models/brm_landreth_fish_watershedpca12345.rds")

brm_landreth_fish_CombinedBiotiapca1 = update(brm_dummy,
                                         newdata = landreth_fish_data,
                                         formula = .~ CombinedBiotiaPC1_s + (1 | stream) + (1 | watershed),
                                         prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
                                                   prior(exponential(5), class = "sd")),
                                         chains = 4, iter = 2000)

saveRDS(brm_landreth_fish_CombinedBiotiapca1, file = "models/brm_landreth_fish_CombinedBiotiapca1.rds")

brm_landreth_fish_CombinedBiotiapca2 = update(brm_dummy,
                                         newdata = landreth_fish_data,
                                         formula = .~ CombinedBiotiaPC2_s + (1 | stream) + (1 | watershed),
                                         prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
                                                   prior(exponential(5), class = "sd")),
                                         chains = 4, iter = 2000)

saveRDS(brm_landreth_fish_CombinedBiotiapca2, file = "models/brm_landreth_fish_CombinedBiotiapca2.rds")

brm_landreth_fish_CombinedBiotiapca12345 = update(brm_dummy,
                                             newdata = landreth_fish_data,
                                             formula = .~ CombinedBiotiaPC1_s + CombinedBiotiaPC2_s + CombinedBiotiaPC3_s +
                                               CombinedBiotiaPC4_s + CombinedBiotiaPC5_s + (1 | stream) + (1 | watershed),
                                             prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
                                                       prior(exponential(5), class = "sd")),
                                             chains = 4, iter = 2000)

saveRDS(brm_landreth_fish_CombinedBiotiapca12345, file = "models/brm_landreth_fish_CombinedBiotiapca12345.rds")
