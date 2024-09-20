library(tidyverse)
library(janitor)
library(tidybayes)

# stream lambdas ----------------------------------------------------------
post_dots_stream = readRDS(file = "posteriors/post_dots.rds") %>%
  filter(grepl("ntercept", model_name)) %>% 
  mutate(taxon_new = case_when(grepl("fishmacros", taxon) ~ "a) fish + macros",
                               grepl("macros", taxon) ~ "b) macros",
                               TRUE ~ "c) fish"))

post_dots_stream %>% 
  group_by(taxon_new, stream) %>% 
  median_qi(.epred) %>% 
  mutate(across(c(.epred, .upper, .lower), ~ round(., 2))) %>% 
  mutate(summary = paste0(.epred, " (", .lower, ", ", .upper, ")")) %>% 
  select(taxon_new, stream, summary) %>% 
  pivot_wider(names_from = taxon_new, values_from = summary) %>% 
  arrange(`a) fish + macros`)

post_dots_stream %>% 
  group_by(taxon_new) %>% 
  median_qi(.epred)

# posterior averaged parameter values ------------------------------------------------------
post_average_parameters = readRDS(file = "posteriors/post_average_parameters_all.rds")

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

post_table_temp = post_average_parameters %>%
  select(-b_Intercept) %>% 
  pivot_longer(cols = starts_with("b_")) %>% 
  mutate(diff = value - 0) %>% 
  group_by(name, taxon) %>% 
  summarize(median = median(value),
            upper = quantile(value, 0.975),
            lower = quantile(value, 0.025),
            prob_greater = sum(diff>0)/max(.draw)) %>% 
  arrange(-prob_greater) %>% 
  left_join(clean_names) 

post_table = post_table_temp %>% 
  select(clean_name, taxon, median, upper, lower) %>% 
  mutate(across(c(median, upper, lower), ~ round(., 2))) %>% 
  mutate(summary = paste0(median, " (", lower, ", ", upper, ")")) %>% 
  select(-median, -upper, -lower, -name) %>% 
  pivot_wider(names_from = taxon, values_from = summary) %>% 
  select(clean_name, `fish + macros`, fish, macros)

prob_table = post_table_temp %>% 
  select(clean_name, taxon, prob_greater, -name) %>% 
  mutate(across(c(prob_greater), ~ round(., 2))) %>% 
  pivot_wider(names_from = taxon, values_from = prob_greater) %>% 
  select(clean_name, `fish + macros`, fish, macros)

write_csv(post_table, file = "tables/post_table.csv")
write_csv(prob_table, file = "tables/prob_table.csv")
