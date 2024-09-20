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
landreth_macros_data = readRDS(file = "data/landreth_macros_data.rds") %>% 
  mutate(counts = as.integer(counts*sample_area_m2))
# check for NAs in trophic groups
# landreth_macros_data %>% group_by(stream) %>% filter(is.na(fishmacros_invertivore_s))

# brm_landreth_macros_topo = brm(dw_g|vreal(counts, xmin, xmax) ~ p_mda_elev_s + p_mda_slope_s + p_da_s +
#                                  (1|stream) + (1|watershed),
#                                data = landreth_macros_data,
#                                stanvars = stanvars,
#                                family = paretocounts(),
#                                prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
#                                          prior(exponential(5), class = "sd"),
#                                          prior(normal(0, 1), class = "b")),
#                                iter = 2000, chains = 4)

saveRDS(brm_landreth_macros_topo, file = "models/brm_landreth_macros_topo.rds")

brm_dummy = readRDS("models/brm_landreth_macros_topo.rds")

brm_elev_std_macros = update(brm_dummy, formula = . ~ p_elva_std_s + (1|stream) + (1|watershed),
                             newdata = landreth_macros_data,
                             iter = 2000, chains = 4,
                             prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
                                       prior(exponential(5), class = "sd")),
                             cores = 4)
saveRDS(brm_elev_std_macros, file = "models/brm_elev_std_macros.rds")

brm_slope_std_macros = update(brm_dummy, formula = . ~ p_slope_std_s + p_mda_slope_s + (1|stream) + (1|watershed),
                              newdata = landreth_macros_data,
                              iter = 2000, chains = 4,
                              prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
                                        prior(exponential(5), class = "sd")),
                              cores = 4)
saveRDS(brm_slope_std_macros, file = "models/brm_slope_std_macros.rds")


brm_herbivore = update(brm_landreth_macros_topo, formula = . ~ fishmacros_herbivore_s + (1|stream) + (1|watershed),
                       newdata = landreth_macros_data,
                       iter = 2000, chains = 4,
                       prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
                                 prior(exponential(5), class = "sd")), 
                       cores = 4)
saveRDS(brm_herbivore, file = "models/brm_macros_herbivore.rds")

brm_invertivore = update(brm_landreth_macros_topo, formula = . ~ fishmacros_invertivore_s + (1|stream) + (1|watershed),
                         newdata = landreth_macros_data,
                         iter = 2000, chains = 4,
                         prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
                                   prior(exponential(5), class = "sd")), 
                         cores = 4)
saveRDS(brm_invertivore, file = "models/brm_macros_invertivore.rds")

brm_omnivore = update(brm_landreth_macros_topo, formula = . ~ fishmacros_omnivore_s + (1|stream) + (1|watershed),
                      newdata = landreth_macros_data,
                      iter = 2000, chains = 4,
                      prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
                                prior(exponential(5), class = "sd")), 
                      cores = 4)
saveRDS(brm_omnivore, file = "models/brm_macros_omnivore.rds")


brm_landreth_macros_intercept = update(brm_landreth_macros_topo,
                                       newdata = landreth_macros_data, formula = . ~ 1 + (1|stream) + (1|watershed),
                                       prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
                                                 prior(exponential(5), class = "sd")))

saveRDS(brm_landreth_macros_intercept, file = "models/brm_landreth_macros_intercept.rds")
