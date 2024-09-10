library(tidybayes)
library(tidyverse)
library(brms)

source("code/custom_functions.R")

# load data and models -----------------------------------------
landreth_fishmacros_data = readRDS(file = "data/landreth_fishmacros_data.rds")
model_list = readRDS(file = "models/temporary/filtered_model_list.rds")

# posterior predictive checks ---------------------------------------------------------

# compare geometric mean
post_gm_list = NULL

for(i in 1:length(model_list)){
  post_gm_list[[i]] = isd_ppcheck(model_list[[i]], re_formula = NULL) %>% 
    mutate(model = as.character(model_list[[i]]$formula)[[1]])
}

post_gm = bind_rows(post_gm_list)

post_gm %>%
  ggplot(aes(x = gm)) +
  geom_histogram(data = . %>% 
                   filter(.draw > 0),
                 bins = 50) +
  scale_x_log10() +
  geom_vline(data = . %>% filter(.draw == 0),
             aes(xintercept = gm), color = 'red') +
  facet_grid(stream_watershed ~ model)

post_gm %>% 
  filter(.draw > 0) %>% 
  ggplot(aes(y = reorder(stream_watershed, median_gm), x = gm)) +
  stat_halfeye() +
  scale_x_log10() +
  geom_point(data = post_gm %>% filter(.draw == 0), color = 'red') +
  facet_wrap(~model)


# bayes_p 
post_gm %>% filter(.draw > 0) %>% 
  left_join(post_gm %>% filter(.draw == 0) %>% 
              rename(y_gm = gm) %>% 
              distinct(stream_watershed, y_gm, model)) %>%
  mutate(diff = gm - y_gm) %>%
  group_by(stream_watershed, model) %>% 
  reframe(bayes_p = sum(diff < 0)/max(.draw)) %>% 
  pivot_wider(names_from = model, values_from = bayes_p)

post_gm %>% filter(.draw > 0) %>% 
  left_join(post_gm %>% filter(.draw == 0) %>% 
              rename(y_gm = gm) %>% 
              distinct(stream_watershed, y_gm, model)) %>%
  mutate(diff = gm - y_gm) %>%
  group_by(stream_watershed, model) %>% 
  reframe(bayes_p = sum(diff < 0)/max(.draw)) %>% 
  group_by(model) %>% 
  median_qi(bayes_p)


