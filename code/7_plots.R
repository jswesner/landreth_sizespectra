library(tidyverse)
library(janitor)
library(tidybayes)
library(brms)
library(ggthemes)
library(isdbayes)
library(viridis)
library(ggh4x)
theme_set(brms::theme_default())

# load data
landreth_fishmacros_data = readRDS(file = "data/landreth_fishmacros_data.rds") %>% mutate(taxon = "fish + macros")
landreth_macros_data = readRDS(file = "data/landreth_macros_data.rds") %>% mutate(taxon = "macros")
landreth_fish_data = readRDS(file = "data/landreth_fish_data.rds") %>% mutate(taxon = "fish")

landreth_data = bind_rows(landreth_fishmacros_data, landreth_macros_data, landreth_fish_data)

# posterior averaged parameter values ------------------------------------------------------
post_average_fishmacros_parameters = readRDS(file = "posteriors/post_average_parameters.rds") %>% mutate(taxon = "fish + macros")
post_average_macros_parameters = readRDS(file = "posteriors/post_average_macros_parameters.rds") %>% mutate(taxon = "macros")
post_average_fish_parameters = readRDS(file = "posteriors/post_average_fish_parameters.rds") %>% mutate(taxon = "fish")

post_average_parameters = bind_rows(post_average_fishmacros_parameters, 
                                    post_average_macros_parameters,
                                    post_average_fish_parameters)

saveRDS(post_average_parameters, file = "posteriors/post_average_parameters_all.rds")


clean_names = post_average_parameters %>% select(starts_with("b_"), taxon) %>% 
  pivot_longer(cols = -taxon) %>% 
  filter(name != "b_Intercept") %>% 
  distinct(name) %>% 
  mutate(clean_name = case_when(grepl("herbivore", name) ~ "% Herbivore Biomass",
                                grepl("omnivore", name) ~ "% Omnivore Biomass",
                                grepl("predator", name) ~ "% Predator Biomass",
                                grepl("p_ag", name) ~ "% Agriculture",
                                grepl("p_slope", name) ~ "Slope SD",
                                grepl("p_da", name) ~ "DA",
                                grepl("p_elva", name) ~ "Elevation SD",
                                grepl("p_watershed", name) ~ "Watershed Size Class"))

posterior_avg_probs = post_average_parameters %>%
  select(-b_Intercept) %>% 
  pivot_longer(cols = starts_with("b_")) %>% 
  mutate(diff = value - 0) %>% 
  group_by(name, taxon) %>% 
  summarize(prob_greater = sum(diff>0)/max(.draw)) %>% 
  arrange(-prob_greater) %>% 
  left_join(clean_names) 

plot_post_avg_parameters = post_average_parameters %>%
  select(-b_Intercept) %>% 
  select(starts_with("b_"), .draw, taxon) %>% 
  pivot_longer(starts_with("b_")) %>% 
  group_by(name) %>% 
  mutate(median_value = median(value),
         dist_0 = abs(median_value) - 0) %>%
  left_join(posterior_avg_probs) %>% 
  mutate(taxon = case_when(taxon == "macros" ~ "b) macros",
                           taxon == "fish" ~ "c) fish",
                           TRUE ~ "a) fish + macros")) %>% 
  ggplot(aes(y = reorder(clean_name, median_value), x = value, fill = median_value)) + 
  stat_slabinterval(color = 'grey50', shape = 1) +
  geom_vline(xintercept = 0) +
  guides(color = "none",
         fill = "none",
         alpha = "none") +
  xlim(-0.1, 0.1) +
  labs(y = "Parameter",
       x = "Posterior Averaged Parameter Value",
       fill = "P(value > 0)") +
  facet_wrap(~taxon, scales = "free_x") +
  NULL

# ggview::ggview(plot_post_avg_parameters, width = 6.5, height = 6.5)
ggsave(plot_post_avg_parameters, width = 6.5, height = 6.5, 
       file = "plots/plot_post_avg_parameters.jpg", dpi = 500)

# posterior averaged regressions ------------------------------------------------------

post_average_parameters = readRDS(file = "posteriors/post_average_parameters_all.rds")
post_dots = readRDS(file = "posteriors/post_dots.rds")

post_medians = post_average_parameters %>%
  select(-b_Intercept) %>% 
  pivot_longer(starts_with("b_")) %>% 
  group_by(name, taxon) %>% 
  summarize(median_effect_size = median(value))

preds = landreth_data %>% select(-dw_g, -counts, -xmin, -xmax,
                                            -sample_area_m2, -n, -watershed, -stream, -taxon) %>% 
  distinct() %>%
  summarise(across(everything(), list(min = min, max = max))) %>%
  pivot_longer(cols = everything()) %>% 
  mutate(min_max = str_sub(name, -3, -1),
         name = str_remove(name, "_max"),
         name = str_remove(name, "_min")) %>%
  pivot_wider(names_from = min_max, values_from = value) %>% 
  rowwise() %>%
  mutate(sequence = list(seq(min, max, length.out = 50))) %>%
  unnest(sequence) %>% 
  select(name, sequence) %>% 
  rename(x = sequence) %>% 
  mutate(name = paste0("b_", name))

post_average_lines = post_average_parameters %>%
  rename(Intercept = b_Intercept) %>%
  pivot_longer(cols = starts_with("b_")) %>% 
  left_join(preds, relationship = "many-to-many") %>% 
  mutate(y = Intercept + x*value) %>% 
  left_join(post_medians) %>% 
  left_join(clean_names)

post_dots_summary = post_dots %>% 
  group_by(value, stream, watershed, name, .draw) %>% 
  summarize(.epred = mean(.epred)) %>% 
  group_by(value, stream, watershed, name) %>% 
  median_qi(.epred) %>% 
  mutate(name = paste0("b_", name)) %>% 
  filter(name != "b_intercept_s") %>% 
  left_join(post_medians) %>% 
  left_join(clean_names)

plot_post_average_lines = post_average_lines %>% 
  filter(.draw <= 500) %>% 
  ggplot(aes(x = x, y = y)) + 
  stat_lineribbon(aes(fill = taxon), .width = c(0.5, 0.75, 0.95), alpha = 0.4) + 
  ggh4x::facet_grid2(reorder(clean_name, -median_effect_size) ~ taxon) +
  ggthemes::scale_fill_colorblind() +
  guides(fill = "none") +
  labs(y = "\u03bb", 
        x= "Predictor (z-score)") +
  theme(strip.text.y = element_text(size = 6))

ggview::ggview(plot_post_average_lines, width = 5, height = 9)
ggsave(plot_post_average_lines, file = "plots/plot_post_average_lines.jpg", width = 5, height = 9, dpi = 300)

# stream lambdas ----------------------------------------------------------

post_dots_stream = readRDS(file = "posteriors/post_dots.rds") %>%
  filter(grepl("ntercept", model_name)) %>% 
  mutate(taxon_new = case_when(grepl("fishmacros", taxon) ~ "a) fish + macros",
                                grepl("macros", taxon) ~ "b) macros",
                                TRUE ~ "c) fish"))

fishmacros_median = post_dots_stream %>% filter(taxon_new == "a) fish + macros") %>% 
  group_by(stream) %>% 
  reframe(median_lambda = median(.epred))

stream_lambdas = post_dots_stream %>%
  left_join(fishmacros_median) %>% 
  ggplot(aes(x = reorder(stream, -median_lambda), y = .epred, color = taxon_new)) + 
  stat_pointinterval(position = position_dodge(width = 0.4),aes(shape = taxon_new),
                     size = 0.01) +
  scale_color_colorblind() +
  labs(x = "Stream Site", 
       y = "\u03bb",
       color = "",
       shape = "") +
  facet_wrap(~taxon_new, ncol = 3) +
  guides(color = "none",
         shape = "none") +
  theme(legend.text = element_text(size = 7),
        text = element_text(size = 9),
        axis.text.x = element_text(angle = 90, hjust = 0))
  
ggsave(stream_lambdas, file = "plots/stream_lambdas.jpg", dpi = 500, width = 6.5, height = 3)


# compare landreth and mle_bin --------------------------------------------
post_summaries_stream = readRDS(file = "posteriors/post_dots.rds") %>%
  filter(grepl("ntercept", model_name)) %>% 
  group_by(stream, taxon) %>% 
  mean_qi(.epred) 

landreth_slopes = read_csv("data/landreth_slopes.csv") %>% 
  pivot_longer(cols = -stream) %>% 
  separate(name, into = c("taxon", "measure")) %>% 
  rename(landreth_mean = value) %>% 
  filter(measure == "slope") %>% 
  mutate(taxon = case_when(taxon == "combined" ~ "fish + macros",
                           TRUE ~ taxon))

post_summaries_stream %>% left_join(landreth_slopes) %>% 
  ggplot(aes(x = landreth_mean, y = .epred)) +
  geom_pointrange(aes(ymin = .lower, ymax = .upper, color = taxon)) +
  facet_wrap(~taxon)


# model comparison plots --------------------------------------------------
fishmacros_waic = readRDS(file = "posteriors/fishmacros_waic.rds")
macros_waic = readRDS(file = "posteriors/macros_waic.rds")
fish_waic = readRDS(file = "posteriors/fish_waic.rds")

all_waic = bind_rows(fishmacros_waic, macros_waic, fish_waic) %>% 
  clean_names() %>% 
  # mutate(model = as.character(formula)) %>%
  mutate(model = str_remove(model_formula, "dw_g \\| vreal\\(counts, xmin, xmax\\) "),
         model = str_replace(model, "\\+ \\(1 \\| stream\\) \\+ \\(1 \\| watershed\\)", "..."),
         model = str_replace(model, " \\(1 \\| stream\\) \\+ \\(1 \\| watershed\\)", " Intercept only"),
         model = as.factor(model),
         model = fct_relevel(model, "~ Intercept only")) %>%
  mutate(lower = waic - 1.96*waic_se,
         upper = waic + 1.96*waic_se) %>% 
  group_by(data) %>% 
  mutate(mean = mean(waic),
         estimate_s = waic/ mean,
         lower_s = lower / mean,
         upper_s = upper / mean) %>% 
  mutate(model_short = case_when(grepl("ntercept", model) ~ "Intercept only",
                                 grepl("invert", model) ~ "Trophic: invertivores",
                                 grepl("omniv", model) ~ "Trophic: omnivores",
                                 grepl("herb", model) ~ "Trophic: herbivores",
                                 grepl("p_devel", model) ~ "Land Use",
                                 TRUE ~ "Topography"))
  

model_selection = all_waic %>% 
  ggplot(aes(y = reorder(model_name, -estimate_s), 
             x = waic, 
             xmin = waic - waic_se, 
             xmax = waic + waic_se)) +
  geom_pointrange(size = 0.1) +
  facet_wrap(~data) +
  labs(x = "WAIC (standardized; lower is better)",
       y = "Model") +
  theme(text = element_text(size = 8))

ggsave(model_selection, file = "plots/model_selection.jpg", width = 6.5, height = 4)
