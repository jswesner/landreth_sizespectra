library(tidyverse)
library(readxl)
library(janitor)

# !!!CHECK IF COUNTS ARE COLLATED PROPERLY IN THE FULL DATA SET (5/21/2024)

# load macros
macro_body_sizes = read_excel("data/body_sizes.xlsx", sheet = "Macros") %>% clean_names() %>% 
  rename(family = famiyl) %>% 
  mutate(sample_id_macros = paste0(site, date),
         sample_id = paste0("macros_", sample_id_macros),
         sample_area_m2 = 0.09*9,
         taxon_group = "macros",
         sample_area_description = "Sample area is a total across 9 replicate surbers. Each surber is 0.09 m2, so total sample area is 0.09*9. From the ms: I used a 0.09 m² Surber Sampler with 500 µm mesh to collect macroinvertebrates from 9 sites that are evenly divided between 1-3 riffles within the reach (depending on the number of riffles available in the reach), according to the NAMC Protocol for the Collection of Aquatic Macroinvertebrate Samples (NAMC 2015).")

# load predictors
abiotic_predictors = read_excel("data/predictors.xlsx") %>% clean_names() %>% 
  rename_with(~paste0("p_", .), everything()) %>% 
  rename(stream = p_stream,
         watershed = p_watershed) %>% 
  select(-contains("p_asin")) %>%  # remove arcsin transformed variables
  pivot_longer(cols = starts_with("p_")) %>% 
  group_by(name) %>% 
  mutate(mean = mean(value, na.rm = T),
         sd = sd(value, na.rm = T),
         value = (value - mean)/sd) %>% 
  mutate(name = paste0(name, "_s")) %>% 
  select(-mean, -sd) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  select(-contains("p_fish_"))

all_percent_trophic_bymass = readRDS(file = "data/all_percent_trophic_bymass.rds") # percent of mass in trophic groups for all, fishonly, and macrosonly
# all_pca_predictors = read_csv(file = "data/all_pca_predictors.csv")

all_predictors = left_join(abiotic_predictors, all_percent_trophic_bymass) %>% 
  # left_join(all_pca_predictors) %>%
  pivot_longer(cols = c(-stream, -watershed, -contains("_s"))) %>% 
  group_by(name) %>% 
  mutate(value = scale(value),
         name = paste0(name, "_s")) %>% 
  pivot_wider(names_from = name, values_from = value) %>%
  select(where(~ !any(is.na(.)))) %>%
  mutate(across(where(~ !is.character(.)), as.numeric)) # removes attributes from scaled variables

write_csv(all_predictors, file = "data/all_predictors.csv")


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

saveRDS(fish_body_sizes, file = "data/fish_body_sizes.rds")

# combine macros and fish and save
landreth_data = bind_rows(fish_body_sizes, macro_body_sizes) %>% 
  group_by(sample_id) %>% 
  mutate(sample_id_no = cur_group_id())

landreth_fish_data = landreth_data %>% 
  filter(taxon_group == "fish") %>%
  group_by(stream) %>% 
  mutate(xmin = min(dw_g),
         xmax = max(dw_g)) %>% 
  group_by(stream, xmin, xmax, dw_g, sample_area_m2) %>% 
  tally(name = "counts") %>% 
  mutate(sample_area_km2 = sample_area_m2/1000) %>% 
  mutate(counts = counts/sample_area_km2) %>% 
  group_by(stream) %>% 
  add_tally() %>% 
  arrange(n) %>% 
  filter(n > 1) %>% 
  ungroup %>% 
  left_join(all_predictors) %>% 
  filter(!is.na(watershed))

landreth_macros_data = landreth_data %>% 
  filter(taxon_group != "fish") %>%
  group_by(stream) %>% 
  mutate(xmin = min(dw_g),
         xmax = max(dw_g)) %>% 
  group_by(stream, xmin, xmax, dw_g, sample_area_m2) %>% 
  tally(name = "counts") %>% 
  mutate(counts = counts/sample_area_m2) %>% 
  group_by(stream) %>% 
  add_tally() %>% 
  arrange(n) %>% 
  filter(n > 1) %>% 
  ungroup %>% 
  left_join(all_predictors) %>% 
  filter(!is.na(watershed))

landreth_fishmacros_data = landreth_data %>% 
  # filter(taxon_group != "fish") %>% 
  group_by(stream) %>% 
  mutate(xmin = min(dw_g),
         xmax = max(dw_g)) %>% 
  group_by(stream, xmin, xmax, dw_g, sample_area_m2) %>% 
  tally(name = "counts") %>% 
  mutate(counts = counts/sample_area_m2) %>% 
  group_by(stream) %>% 
  add_tally() %>% 
  arrange(n) %>% 
  filter(n > 1) %>% 
  ungroup %>% 
  left_join(all_predictors) %>% 
  filter(!is.na(watershed))

saveRDS(landreth_fish_data, file = "data/landreth_fish_data_uncorrected.rds")
saveRDS(landreth_macros_data, file = "data/landreth_macros_data_uncorrected.rds")
saveRDS(landreth_fishmacros_data, file = "data/landreth_fishmacros_data_uncorrected.rds")

peak_min_sizes_macros <- readRDS("data/peak_min_sizes_macros.RDS")
peak_min_sizes_fish = readRDS("data/peak_min_sizes_fish.rds")
peak_min_sizes_fishmacros = readRDS("data/peak_min_sizes.rds")

landreth_macros_data = readRDS(file = "data/landreth_macros_data_uncorrected.rds") %>% 
  left_join(peak_min_sizes_macros) %>% 
  filter(dw_g >= sum_min) %>% 
  group_by(stream) %>% 
  mutate(xmin = min(dw_g),
         xmax = max(dw_g))

landreth_fishmacros_data = readRDS(file = "data/landreth_fishmacros_data_uncorrected.rds") %>% 
  left_join(peak_min_sizes_fishmacros) %>% 
  filter(dw_g >= sum_min) %>% 
  group_by(stream) %>% 
  mutate(xmin = min(dw_g),
         xmax = max(dw_g))

landreth_fish_data = readRDS(file = "data/landreth_fish_data_uncorrected.rds") %>% 
  left_join(peak_min_sizes_fish) %>% 
  filter(dw_g >= sum_min) %>% 
  group_by(stream) %>% 
  mutate(xmin = min(dw_g),
         xmax = max(dw_g))

saveRDS(landreth_macros_data, file = "data/landreth_macros_data.rds")
saveRDS(landreth_fish_data, file = "data/landreth_fish_data.rds")
saveRDS(landreth_fishmacros_data, file = "data/landreth_fishmacros_data.rds")




# check correlations
library(ggcorrplot)

correlations = all_predictors %>% 
  select(c(-stream, -watershed)) %>%
  select(starts_with("Watershed")) %>% 
  cor()

ggcorrplot(correlations)

landreth_data %>% glimpse %>% 
  group_by(stream, date) %>% 
  reframe(mean = mean(dw_g))

all_predictors %>% 
  select(starts_with("fishmacros")) %>% 
  select(!ends_with("_s")) %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(x = value, fill = name)) + 
  geom_histogram() +
  facet_wrap(~name)

               