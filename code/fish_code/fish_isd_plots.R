library(tidyverse)
library(janitor)
library(tidybayes)
library(brms)
library(ggthemes)
library(isdbayes)
library(viridis)
theme_set(brms::theme_default())

# plots will save here
directory = "plots/fish_plots/"

#1) Load models and data
# load data and models -----------------------------------------
model_files_fish = list.files(path = "models/fish_models", pattern = "\\.rds$", full.names = T)
model_lists_fish <- list()

for (file in model_files_fish) {
  model_list_fish[[file]] <- readRDS(file)
}

dat = as_tibble(model_list_fish[[2]]$data)

#2) Resample data
n_samples = 5000

dat_resampled_rank = dat %>% 
  group_by(stream, watershed) %>% 
  sample_n(n_samples, weight = counts, replace = T) %>% 
  select(stream, counts, watershed, ends_with("_s"), dw_g) %>% 
  group_by(stream) %>% 
  mutate(xmin = min(dw_g),
         xmax = max(dw_g),
         data = "y_raw") %>% 
  arrange(-dw_g) %>% 
  mutate(n_yx = 1:max(row_number())) %>% 
  ungroup()


#3) ISD by sample ----------------------------------

epred_draws = dat_resampled_rank %>% 
  select(stream, watershed, xmin, xmax, ends_with("_s")) %>% 
  distinct() %>% 
  mutate(counts = 1) %>% 
  add_epred_draws(model_list_fish[[2]], re_formula = NULL)

epred_draws_summary = epred_draws %>% 
  group_by(stream, watershed, xmin, xmax) %>% 
  median_qi(.epred) %>% 
  pivot_longer(cols = c(.epred, .lower, .upper),
               names_to = "quantile", 
               values_to = "lambda")

epred_draws_list = epred_draws_summary %>% 
  group_by(stream, quantile, xmin, xmax) %>% 
  group_split()

temp = lapply(epred_draws_list, function(df) {
  df %>% expand_grid(x = 10^seq(log10(xmin), log10(xmax), length.out = 50)) %>% 
    mutate(prob_yx = (1 - (x^(lambda + 1) - (xmin^(lambda + 1))) / ((xmax)^(lambda + 1) - (xmin^(lambda + 1)))),
           n_yx = prob_yx * n_samples) 
})


site_ordered = dat_resampled_rank %>% 
  ungroup %>% 
  distinct(stream) %>%
  arrange(stream) %>% 
  mutate(site_order = 1:nrow(.),
         site_order_site = paste0(site_order, ") ", stream))

isd_lines = bind_rows(temp) %>% 
  left_join(site_ordered)

sample_max = length(unique(dat$stream))

dat_sliced = dat_resampled_rank %>% 
  # filter(sample_id <= sample_max) %>%
  group_by(stream) %>% 
  slice_sample(n = 1000) %>% 
  left_join(site_ordered)

plot_isd = isd_lines %>% 
  # filter(sample_id <= sample_max) %>%
  select(-prob_yx, -lambda) %>% 
  pivot_wider(names_from = quantile, values_from = n_yx) %>% 
  mutate(dw_g = x) %>%
  ggplot(aes(x = dw_g, y = .epred)) +
  geom_point(data = dat_sliced ,
             aes(y = n_yx,
                 color = stream),
             size = 0.5,
             shape = 1)  +
  geom_line() +
  # geom_text(data = site_ordered, aes(label = site_order_site),
  #           x = -1.1,
  #           y = 0.5,
  #           size = 1.2) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = 0.4) +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(ylim = c(0.5, NA)) +  
  facet_wrap(~stream) +
  guides(size = "none",
         color = "none") +
  # theme(axis.text = element_blank(),
  #       strip.text = element_blank(),
  #       panel.spacing = unit(0.1,'lines')) +
  labs(y = "\u03bb",
       x = "gDM") +
  NULL

ggview::ggview(plot_isd, width = 6.5, height = 7)
ggsave(plot_isd, width = 6.5, height = 7,
       file = paste0(directory, "fish_plot_isd.jpg"), dpi = 300)
saveRDS(plot_isd, file = paste0(directory, "fish_plot_isd.rds"))


