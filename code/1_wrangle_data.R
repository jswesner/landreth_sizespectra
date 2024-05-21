library(tidyverse)
library(readxl)
library(janitor)

# load macros
macro_body_sizes = read_excel("data/body_sizes.xlsx", sheet = "Macros") %>% clean_names() %>% 
  rename(family = famiyl) %>% 
  mutate(sample_id_macros = paste0(site, date),
         sample_id = paste0("macros_", sample_id_macros),
         sample_area_m2 = 0.09*9,
         taxon_group = "macros",
         sample_area_description = "Sample area is a total across 9 replicate surbers. Each surber is 0.09 m2, so total sample area is 0.09*9. From the ms: I used a 0.09 m² Surber Sampler with 500 µm mesh to collect macroinvertebrates from 9 sites that are evenly divided between 1-3 riffles within the reach (depending on the number of riffles available in the reach), according to the NAMC Protocol for the Collection of Aquatic Macroinvertebrate Samples (NAMC 2015).")

# load predictors
predictors = read_excel("data/predictors.xlsx") %>% clean_names() %>% 
  rename_with(~paste0("p_", .), everything()) %>% 
  rename(stream = p_stream,
         watershed = p_watershed) %>% 
  pivot_longer(cols = starts_with("p_")) %>% 
  group_by(name) %>% 
  mutate(mean = mean(value, na.rm = T),
         sd = sd(value, na.rm = T),
         value = (value - mean)/sd) %>% 
  mutate(name = paste0(name, "_s")) %>% 
  select(-mean, -sd) %>% 
  pivot_wider(names_from = name, values_from = value)

# load fish
electrofishing_info = "data/electrofishing_info.xlsx" %>% 
  excel_sheets() %>% 
  set_names() %>% 
  map(read_excel, path = "data/electrofishing_info.xlsx")

fish_sample_area = electrofishing_info$stream_info %>% 
  clean_names() %>% 
  select(site, date, length_reach, starts_with("wetted")) %>% 
  mutate(mean_wetted_width_m = (wetted_width1_m + wetted_width2 + wetted_width3)/3) %>% 
  mutate(sample_area_m2 = length_reach*mean_wetted_width_m) %>% 
  select(site, sample_area_m2)

fish_body_sizes = read_excel("data/body_sizes.xlsx", sheet = "Fish") %>% clean_names() %>% 
  group_by(site) %>% 
  mutate(sample_id_fish = paste0(site, date),
         sample_id = paste0("fish_", sample_id_fish),
         taxon_group = "fish",
         dw_mg = dw_g*1000) %>% 
  left_join(fish_sample_area)

# combine macros and fish and save
landreth_data = left_join(bind_rows(fish_body_sizes, macro_body_sizes), predictors) %>% 
  group_by(sample_id) %>% 
  mutate(sample_id_no = cur_group_id())

write_csv(landreth_data, file = "data/landreth_data.csv")

landreth_fish_data = landreth_data %>% 
  filter(taxon_group == "fish") %>% 
  group_by(sample_id) %>% 
  mutate(xmin = min(dw_g),
         xmax = max(dw_g)) %>% 
  group_by(sample_id, xmin, xmax, dw_g, sample_id_no, site, across(starts_with("p_"))) %>% 
  tally(name = "counts") %>% 
  group_by(sample_id) %>% 
  add_tally() %>% 
  arrange(n) %>% 
  filter(n > 1) %>%  # Removes 1 row of a data for a sample_id with just 1 fish. All others have >= 16 fish.
  ungroup

landreth_macros_data = landreth_data %>% 
  filter(taxon_group != "fish") %>% 
  group_by(sample_id) %>% 
  mutate(xmin = min(dw_g),
         xmax = max(dw_g)) %>% 
  group_by(sample_id, xmin, xmax, dw_g, sample_id_no, across(starts_with("p_"))) %>% 
  tally(name = "counts") %>% 
  group_by(sample_id) %>% 
  add_tally() %>% 
  arrange(n) %>% 
  filter(n > 1) %>% 
  ungroup

landreth_fishmacros_data = landreth_data %>% 
  # filter(taxon_group != "fish") %>% 
  group_by(sample_id) %>% 
  mutate(xmin = min(dw_g),
         xmax = max(dw_g)) %>% 
  group_by(sample_id, xmin, xmax, dw_g, sample_id_no, sample_area_m2, across(starts_with("p_"))) %>% 
  tally(name = "counts") %>% 
  group_by(sample_id) %>% 
  add_tally() %>% 
  arrange(n) %>% 
  filter(n > 1) %>% 
  ungroup 

write_csv(landreth_fish_data, file = "data/landreth_fish_data.csv")
write_csv(landreth_macros_data, file = "data/landreth_macros_data.csv")
write_csv(landreth_fishmacros_data, file = "data/landreth_fishmacros_data.csv")


# check correlations
library(ggcorrplot)

correlations = predictors %>% 
  select(starts_with("p_")) %>%
  cor()

ggcorrplot(correlations)

plot(p_elva_std_s ~ p_mda_elev_s, data = predictors)
