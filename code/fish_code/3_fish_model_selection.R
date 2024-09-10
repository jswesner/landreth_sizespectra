library(tidybayes)
library(tidyverse)
library(brms)
library(isdbayes)

source("code/custom_functions.R")

# load data and models -----------------------------------------
landreth_fish_data = readRDS(file = "data/landreth_fish_data.rds")

model_fish_files = list.files(
  path = "models/fish_models/", 
  pattern = "\\.rds$", 
  full.names = TRUE
)

# Exclude "model_fish_list.rds". Otherwise the code below gives an error
model_fish_files = model_fish_files[basename(model_fish_files) != "model_fish_list.rds"]

model_fish_list = NULL

# Loop through each file and read the RDS into the list
for (file in model_fish_files) {
  model_fish_list[[file]] <- readRDS(file)
}


# get waic for intercept only

for(i in 1:length(model_fish_list)){
  model_fish_waic = get_waic_isd(model_fish_list[[i]], nsamples = 1000, ndraws = 500, resp = "dw_g") %>% 
    mutate(model_fish_name = names(model_fish_list[i]),
           model_fish_formula = as.character(model_fish_list[[i]]$formula)[1])
  
  model_fish_list[[i]]$criteria = model_fish_waic
}

saveRDS(model_fish_list, file = "posteriors/model_fish_list.rds")


# posterior averaging -----------------------------------------------------
model_fish_list = readRDS(file = "posteriors/model_fish_list.rds")

# filter to only subset of models
# filtered_model_fish_names <- names(model_fish_list)[!grepl("PC", names(model_fish_list), ignore.case = TRUE)]
filtered_model_fish_names <- names(model_fish_list)[grepl("agtopo|trophic|intercept", names(model_fish_list), ignore.case = TRUE)]
filtered_model_fish_list <- model_fish_list[filtered_model_fish_names]
saveRDS(filtered_model_fish_list, file = "models/fish_models/temporary/filtered_model_fish_list.rds")

# get model weights (adpated from https://github.com/paul-buerkner/brms/blob/master/R/model_weights.R)
waic_values = numeric(length(filtered_model_fish_list)) 

for(i in 1:length(filtered_model_fish_list)) waic_values[i] = filtered_model_fish_list[[i]]$criteria$waic

se_values = numeric(length(filtered_model_fish_list))
for(i in 1:length(filtered_model_fish_list)) se_values[i] = filtered_model_fish_list[[i]]$criteria$waic_se
min_waic <- min(waic_values)
adjusted_diffs <- (waic_values - min_waic) / se_values
weights <- exp(-adjusted_diffs / 2)
model_fish_weights <- weights / sum(weights)  # Normalize to sum to 1

waic_weights = tibble(waic = waic_values, model_fish_weights = model_fish_weights, model = names(filtered_model_fish_list))

# add model weights to data2 in models
for(i in 1:length(filtered_model_fish_list)) filtered_model_fish_list[[i]]$data2 = list(model_fish_weight = model_fish_weights[i])

model_fish_posts = list()

for(i in 1:length(filtered_model_fish_list)){
  model_fish_posts[[i]] = as_draws_df(filtered_model_fish_list[[i]]) %>% 
    mutate(model = names(filtered_model_fish_list[i]),
           weight = filtered_model_fish_list[[i]]$data2$model_fish_weight)
}

post_average_fish_parameters = bind_rows(model_fish_posts) %>% 
  # filter(.draw <= 1000) %>% 
  select(.draw, starts_with("b_"), model, weight) %>% 
  mutate(across(starts_with("b_"), ~ replace_na(.x, 0))) %>%
  group_by(.draw) %>% 
  summarize(across(starts_with("b_"), ~ weighted.mean(.x, weight))) %>% 
  mutate(data = "fish + fish")

saveRDS(post_average_fish_parameters, file = "posteriors/post_average_fish_parameters.rds")


# plot 
post_average_fish_parameters %>%
  select(-b_Intercept) %>% 
  pivot_longer(starts_with("b_")) %>% 
  group_by(name) %>% 
  mutate(median_value = median(value), na.rm = T) %>% 
  ggplot(aes(y = reorder(name, median_value), x = value, color = median_value)) + 
  stat_pointinterval() +
  geom_vline(xintercept = 0)

