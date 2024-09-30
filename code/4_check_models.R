library(tidybayes)
library(tidyverse)
library(brms) 
theme_set(brms::theme_default())

source("code/custom_functions.R")

# load data and models -----------------------------------------
# landreth_fishmacros_data = readRDS(file = "data/landreth_fish_data.rds")
model_list = readRDS(file = "posteriors/model_list.rds")
dat = model_list[[1]]$data # body sizes are the same regardless of model index 1, 2, or 3, but predictors will differ


# 1) resample data
dat_sampled = dat %>% 
  left_join(model_list[[3]]$data %>% distinct(stream, fishmacros_herbivore_s, fishmacros_omnivore_s, fishmacros_predator_s)) %>% 
  sample_n(2000, weight = counts, replace = T)

# 2) get posterior predictions
post_samples_list = list()
for(i in 1:length(model_list)){
  post_samples_list[[i]] = dat_sampled %>% add_predicted_draws(model_list[[i]], ndraws = 100, 
                                                               re_formula = NULL) %>% 
    mutate(model = names(model_list[i]),
           model = str_remove(model, "models/brm_landreth_macros_"))
}
# saveRDS(post_samples_list, file = "posteriors/post_samples_list.rds")
# post_samples_list = readRDS(file = "posteriors/post_samples_list.rds")
post_samples = bind_rows(post_samples_list) 
dat_expanded = dat_sampled %>% expand_grid(model = unique(post_samples$model)) # make model-specific data

# 3) pp_check densities
post_samples %>% 
  filter(.draw <= 20) %>% 
  ggplot(aes(x = .prediction)) +
  geom_density(aes(group = .draw)) +
  scale_x_log10() +
  facet_wrap(~model) + 
  geom_density(data = dat_expanded, aes(x = dw_g), 
               color = "dodgerblue")

# 4) pp_check stats 
post_median_gm = post_samples %>% 
  group_by(.draw, model) %>% 
  reframe(median = median(.prediction),
          gm = exp(mean(log(.prediction)))) %>% 
  pivot_longer(cols = c(median, gm))

dat_median_gm = dat_expanded %>% 
  group_by(model) %>% 
  reframe(median = median(dw_g),
          gm = exp(mean(log(dw_g)))) %>% 
  pivot_longer(cols = c(median, gm))

post_median_gm %>% 
  ggplot(aes(x = value)) + 
  geom_histogram(bins = 50) +
  scale_x_log10() +
  ggh4x::facet_grid2(model ~ name) +
  geom_vline(data = dat_median_gm, aes(xintercept = value))


# 5) bayesian p-values

post_median_gm %>% 
  left_join(dat_median_gm %>% rename(value_raw = value)) %>% 
  mutate(diff = value_raw - value) %>% 
  group_by(name, model) %>% 
  reframe(bayes_p = sum(diff>0)/max(.draw))



