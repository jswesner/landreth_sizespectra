library(tidyverse)
library(brms)
library(isdbayes)

# model fish only ---------------------------------------------------------
# load data
landreth_fish_data = read_csv(file = "data/landreth_fish_data.csv") 

# plot
landreth_fish_data %>% 
  uncount(weights = counts) %>% 
  group_by(sample_id) %>% 
  arrange(sample_id, -dw_g) %>% 
  mutate(order = row_number()) %>% 
  ggplot(aes(x = dw_g, y = order)) + 
  geom_point(shape = 1) +
  scale_x_log10() + 
  scale_y_log10()

brm_landreth_fish = brm(dw_g|vreal(counts, xmin, xmax) ~ 1 + (1|sample_id),
                        data = landreth_fish_data,
                        stanvars = stanvars,
                        family = paretocounts(),
                        prior = c(prior(normal(-1.2, 0.8), class = "Intercept"),
                                  prior(exponential(5), class = "sd")),
                        iter = 2000, chains = 4, cores = 4)

saveRDS(brm_landreth_fish, file = "models/brm_landreth_fish.rds")

brm_landreth_fish = readRDS(file = "models/brm_landreth_fish.rds")
