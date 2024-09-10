library(tidyverse)
library(readxl)
library(janitor)

landreth_macros_data_ffg <- readxl::read_excel("data/landreth_fishmacros_data_ffg.xlsx", 
                                           sheet = "BMI FFG") %>% 
  mutate(trophic = case_when(FFG_Score == 1 ~ "herbivore",
                             FFG_Score == 2 ~ "omnivore", 
                             FFG_Score == 3 ~ "predator")) %>% 
  rename(family = Family) %>% 
  clean_names()

fish_body_sizes = readRDS(file = "data/fish_body_sizes.rds")

# In the original landreth fish ffg file, rock bass had two entries (predator and invertivore)
# I deleted the predator entry before loading here, since fish cannot have more than one trophic assignment.
# I also added trophic assignments for potential hybrids that were the same as the first hybrid species.
landreth_fish_ffg = read_excel("data/Landreth Fish FFGs.xlsx", sheet = "FFG_work") %>% 
  clean_names() 

fish_trophic = read_csv("data/fish_trophic.csv") %>% 
  arrange(scientific) %>% 
  select(-trophic_gpt, -trophic_fishbase, -trophic_wesner) %>% 
  left_join(landreth_fish_ffg) %>% 
  mutate(trophic = case_when(grepl("invert", ffg) ~ "invertivore", TRUE ~ ffg),
         trophic = tolower(trophic))


fish_trophic_size = fish_body_sizes %>% 
  left_join(fish_trophic) %>% 
  select(stream, dw_g, sample_id, sample_area_m2, trophic, ffg) %>% 
  group_by(site, stream, sample_id, trophic, ffg, sample_area_m2) %>% 
  mutate(dw_g_m2 = dw_g/sample_area_m2) %>% 
  mutate(taxon = "fish")

macro_body_sizes = read_excel("data/body_sizes.xlsx", sheet = "Macros") %>% clean_names() %>% 
  rename(family = famiyl) %>% 
  mutate(sample_id_macros = paste0(site, date),
         sample_id = paste0("macros_", sample_id_macros),
         sample_area_m2 = 0.09*9,
         taxon_group = "macros",
         sample_area_description = "Sample area is a total across 9 replicate surbers. Each surber is 0.09 m2, so total sample area is 0.09*9. From the ms: I used a 0.09 m² Surber Sampler with 500 µm mesh to collect macroinvertebrates from 9 sites that are evenly divided between 1-3 riffles within the reach (depending on the number of riffles available in the reach), according to the NAMC Protocol for the Collection of Aquatic Macroinvertebrate Samples (NAMC 2015).")


macro_trophic_size = macro_body_sizes %>% 
  left_join(landreth_macros_data_ffg) %>% 
  select(site, stream, dw_g, sample_id, sample_area_m2, trophic, ffg) %>% 
  mutate(dw_g_m2 = dw_g/sample_area_m2) %>% 
  mutate(taxon = "macros")

stream_trophic = bind_rows(fish_trophic_size, macro_trophic_size) %>% 
  ungroup %>% 
  select(stream, trophic, dw_g_m2) %>% 
  group_by(stream) %>% 
  mutate(total_dw_g_m2 = sum(dw_g_m2)) %>% 
  group_by(stream, trophic, total_dw_g_m2) %>% 
  reframe(sum_dw_g_m2_trophic = sum(dw_g_m2)) %>% 
  mutate(perc_trophic = sum_dw_g_m2_trophic/total_dw_g_m2) %>% 
  select(stream, trophic, perc_trophic) %>% 
  mutate(trophic = paste0("fishmacros_", trophic)) %>% 
  pivot_wider(names_from = trophic, values_from = perc_trophic)

fish_trophic = bind_rows(fish_trophic_size, macro_trophic_size) %>% 
  ungroup %>% 
  filter(taxon == "fish") %>% 
  select(stream, trophic, dw_g_m2) %>% 
  group_by(stream) %>% 
  mutate(total_dw_g_m2 = sum(dw_g_m2)) %>% 
  group_by(stream, trophic, total_dw_g_m2) %>% 
  reframe(sum_dw_g_m2_trophic = sum(dw_g_m2)) %>% 
  mutate(perc_trophic = sum_dw_g_m2_trophic/total_dw_g_m2) %>% 
  select(stream, trophic, perc_trophic) %>% 
  mutate(trophic = paste0("fishonly_", trophic)) %>% 
  pivot_wider(names_from = trophic, values_from = perc_trophic)

macro_trophic = bind_rows(fish_trophic_size, macro_trophic_size) %>% 
  ungroup %>% 
  filter(taxon != "fish") %>% 
  select(stream, trophic, dw_g_m2) %>% 
  group_by(stream) %>% 
  mutate(total_dw_g_m2 = sum(dw_g_m2)) %>% 
  group_by(stream, trophic, total_dw_g_m2) %>% 
  reframe(sum_dw_g_m2_trophic = sum(dw_g_m2)) %>% 
  mutate(perc_trophic = sum_dw_g_m2_trophic/total_dw_g_m2) %>% 
  select(stream, trophic, perc_trophic) %>% 
  mutate(trophic = paste0("macrosonly_", trophic)) %>% 
  pivot_wider(names_from = trophic, values_from = perc_trophic)

all_percent_trophic_bymass = left_join(stream_trophic, macro_trophic) %>% 
  left_join(fish_trophic)

saveRDS(all_percent_trophic_bymass, file = "data/all_percent_trophic_bymass.rds")
