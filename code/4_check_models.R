library(tidybayes)
library(tidyverse)
library(brms) 
theme_set(brms::theme_default())

source("code/custom_functions.R")

# load data and models -----------------------------------------
# landreth_fishmacros_data = readRDS(file = "data/landreth_fish_data.rds")
model_list = readRDS(file = "posteriors/model_list.rds")
dat = model_list[[1]]$data # body sizes are the same regardless of model index 1, 2, or 3, but predictors will differ


# prior predictive --------------------------------------------------------

prior_post_all <- map_dfr(1:3, function(i) {
  posts <- as_draws_df(model_list[[i]]) %>%
    select(starts_with("sd_"), starts_with("b_"), .draw, "Intercept") %>%
    mutate(prior_post = "posterior") %>%
    select(-b_Intercept)
  
  priors <- tibble(cols = names(posts)) %>%
    mutate(mean = case_when(
      cols == "Intercept" ~ -1.6,
      grepl("^b_", cols) ~ 0,
      grepl("^sd_", cols) ~ 5
    )) %>%
    filter(!is.na(mean)) %>%
    mutate(sd = case_when(
      cols == "Intercept" ~ 0.5,
      grepl("^b_", cols) ~ 0.2
    )) %>%
    expand_grid(.draw = 1:1000) %>%
    mutate(prior_sample = case_when(
      grepl("^sd", cols) ~ rexp(nrow(.), mean),
      TRUE ~ rnorm(nrow(.), mean, sd)
    )) %>%
    select(cols, prior_sample, .draw) %>%
    pivot_wider(names_from = cols, values_from = prior_sample) %>%
    mutate(prior_post = "prior")
  
  bind_rows(posts, priors) %>% mutate(model = paste0("model_", i)) %>% 
    pivot_longer(cols = c(-.draw, -prior_post, -model)) %>% 
    mutate(name == as.factor(name),
           name = fct_relevel(name, 'Intercept'))
})


prior_post_intercept = prior_post_all %>% 
  filter(model == "model_2")  %>% 
  filter(!is.na(value)) %>% 
  ggplot(aes(x = value, fill = prior_post, alpha = prior_post)) +
  geom_density() +
  scale_alpha_manual(values = c(1, 0.4)) +
  ggh4x::facet_wrap2(~name, ncol = 1) +
  ggthemes::scale_fill_colorblind() +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(hjust = 0,
                                  size = 8),
        legend.position = "top") +
  labs(y = "",
       x = "Parameter Value",
       fill = "",
       alpha = "",
       subtitle = "a) Intercept Model")


prior_post_landuse = prior_post_all %>% 
  filter(model == "model_1")  %>% 
  filter(!is.na(value)) %>% 
  ggplot(aes(x = value, fill = prior_post, alpha = prior_post)) +
  geom_density() +
  scale_alpha_manual(values = c(1, 0.4)) +
  ggh4x::facet_wrap2(~name, ncol = 1) +
  ggthemes::scale_fill_colorblind() +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(hjust = 0,
                                  size = 8)) +
  labs(y = "",
       x = "Parameter Value",
       fill = "",
       alpha = "",
       subtitle = "b) Land Use Model") +
  guides(fill = "none", 
         alpha = "none")

prior_post_trophic = prior_post_all %>% 
  filter(model == "model_3")  %>% 
  filter(!is.na(value)) %>% 
  ggplot(aes(x = value, fill = prior_post, alpha = prior_post)) +
  geom_density() +
  scale_alpha_manual(values = c(1, 0.4)) +
  ggh4x::facet_wrap2(~name, ncol = 1) +
  ggthemes::scale_fill_colorblind() +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(hjust = 0,
                                  size = 8),
        legend.position = "top") +
  labs(y = "",
       x = "Parameter Value",
       fill = "",
       alpha = "",
       subtitle = "c) Trophic Model") +
  guides(fill = "none", 
         alpha = "none")


library(patchwork)

prior_post_plot = prior_post_intercept + prior_post_landuse + prior_post_trophic

ggsave(prior_post_plot, file = "plots/prior_post_plot.jpg", width = 6.5, height = 8,
       dpi = 400)



# posterior predictive ----------------------------------------------------


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
post_gm = post_samples %>% 
  group_by(.draw, model) %>% 
  reframe(gm = exp(mean(log(.prediction)))) %>% 
  mutate(model = case_when(grepl("agtopo", model) ~ "b) Land Use Model",
                           grepl("trophic", model) ~ "c) Trophic Model",
                           TRUE ~ "a) Intercept Only Model"))

dat_gm = dat_expanded %>% 
  group_by(model) %>% 
  reframe(gm = exp(mean(log(dw_g)))) %>% 
  mutate(model = case_when(grepl("agtopo", model) ~ "b) Land Use Model",
                           grepl("trophic", model) ~ "c) Trophic Model",
                           TRUE ~ "a) Intercept Only Model"))

labels = tibble(model = "a) Intercept Only Model") %>% 
  mutate(label_1 = "Simulations above the\nempirical value",
         label_2 = "Simulations below the\nempirical value") %>% 
  pivot_longer(cols = starts_with("label")) %>% 
  mutate(gm = c(0.00090, 0.00079),
         hist_color = c("above", "below"))

posterior_predictions = post_gm %>% 
  left_join(dat_gm %>% rename(gm_empirical = gm)) %>% 
  mutate(hist_color = case_when(gm >= gm_empirical ~ "above",
                                TRUE ~ "below")) %>% 
  ggplot(aes(x = gm)) + 
  geom_histogram(bins = 50, aes(fill = hist_color)) +
  scale_x_log10(limits = c(0.00075, 0.001)) +
  ggh4x::facet_wrap2(~model, ncol = 1) +
  geom_vline(data = dat_gm, aes(xintercept = gm)) +
  labs(x = "Geometric mean individual body size (g dry mass)") +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(hjust = 0)) +
  guides(fill = "none",
         color = "none") +
  scale_fill_brewer(type = "qual") +
  scale_color_brewer(type = "qual") +
  geom_text(data = labels, aes(label = value, y = 8, color = hist_color),
            size = 3, hjust = 0) 

ggsave(posterior_predictions, file = "plots/posterior_predictions.jpg", width = 5, height = 8, dpi = 400)


# 5) bayesian p-values

post_gm %>% 
  left_join(dat_gm %>% rename(value_raw = gm)) %>% 
  mutate(diff = value_raw - gm) %>% 
  group_by(model) %>% 
  reframe(bayes_p = sum(diff>0)/max(.draw))



