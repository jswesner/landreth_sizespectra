library(tidybayes)
library(tidyverse)
library(brms)


# fish only ---------------------------------------------------------------
landreth_fish_data = read_csv(file = "data/landreth_fish_data.csv") 
brm_landreth_fish = readRDS("models/brm_landreth_fish.rds")


sample_lambdas = landreth_fish_data %>% 
  distinct(xmin, xmax, sample_id, across(starts_with("p_"))) %>% 
  mutate(counts = 1) %>% 
  add_epred_draws(brm_landreth_fish)




sample_lambdas %>% 
  group_by(sample_id) %>% 
  mutate(median = median(.epred)) %>% 
  ggplot(aes(x = .epred, y = reorder(sample_id, median))) + 
  stat_halfeye()


sample_lambdas %>% 
  group_by(sample_id) %>% 
  mutate(median = median(.epred)) %>% 
  ggplot(aes(x = p_fish_piscivore, y = .epred)) + 
  stat_pointinterval(aes(group = sample_id)) +
  scale_x_log10() +
  NULL


# fish and macros ---------------------------------------------------------
landreth_fishmacros_data = read_csv(file = "data/landreth_fishmacros_data.csv")
brm_landreth_fishmacros = readRDS("models/brm_landreth_fishmacros.rds")

brm_landreth_fishmacros

post_lambdas_fishmacros = brm_landreth_fishmacros$data %>% 
  glimpse %>% 
  distinct(sample_id, xmin, xmax) %>% 
  mutate(counts = 1) %>% 
  add_epred_draws(brm_landreth_fishmacros) %>% 
  left_join(landreth_fishmacros_data %>% select(-xmin, -xmax, -dw_g, -sample_area_m2) %>% distinct())

post_lambdas_fishmacros %>% 
  ggplot(aes(y = sample_id, x = .epred)) + 
  stat_halfeye() 

post_lambdas_fishmacros %>% glimpse %>% 
  ggplot(aes(x = p_forest, y = .epred)) + 
  stat_halfeye(aes(group = sample_id)) +
  geom_smooth(method = lm)
  
