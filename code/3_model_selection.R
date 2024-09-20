library(tidybayes)
library(tidyverse)
library(brms)
library(isdbayes)

source("code/custom_functions.R")

# load data and models -----------------------------------------
landreth_fishmacros_data = readRDS(file = "data/landreth_fishmacros_data.rds")

model_files = list.files(
  path = "models/", 
  pattern = "\\.rds$", 
  full.names = TRUE
)

# Exclude "model_list.rds". Otherwise the code below gives an error
model_files = model_files[basename(model_files) != "model_list.rds"]

model_list = NULL

# Loop through each file and read the RDS into the list
for (file in model_files) {
  model_list[[file]] <- readRDS(file)
}


# get waic for intercept only

for(i in 1:length(model_list)){
  model_waic = get_waic_isd(model_list[[i]], nsamples = 1000, ndraws = 500, resp = "dw_g") %>% 
    mutate(model_name = names(model_list[i]),
           model_formula = as.character(model_list[[i]]$formula)[1])
  
  model_list[[i]]$criteria = model_waic
}

saveRDS(model_list, file = "posteriors/model_list.rds")


# posterior averaging -----------------------------------------------------
model_list = readRDS(file = "posteriors/model_list.rds")

# filter to only subset of models
# filtered_model_names <- names(model_list)[!grepl("PC", names(model_list), ignore.case = TRUE)]
filtered_model_names <- names(model_list)[grepl("agtopo|trophic|intercept", names(model_list), ignore.case = TRUE)]
filtered_model_list <- model_list[filtered_model_names]
saveRDS(filtered_model_list, file = "models/temporary/filtered_model_list.rds")

# get model weights (adpated from https://github.com/paul-buerkner/brms/blob/master/R/model_weights.R)
waic_values = numeric(length(filtered_model_list)) 

for(i in 1:length(filtered_model_list)) waic_values[i] = filtered_model_list[[i]]$criteria$waic

se_values = numeric(length(filtered_model_list))
for(i in 1:length(filtered_model_list)) se_values[i] = filtered_model_list[[i]]$criteria$waic_se
min_waic <- min(waic_values)
adjusted_diffs <- (waic_values - min_waic) / se_values
weights <- exp(-adjusted_diffs / 2)
model_weights <- weights / sum(weights)  # Normalize to sum to 1

waic_weights = tibble(waic = waic_values, model_weights = model_weights, model = names(filtered_model_list))

# add model weights to data2 in models
for(i in 1:length(filtered_model_list)) filtered_model_list[[i]]$data2 = list(model_weight = model_weights[i])

model_posts = list()

for(i in 1:length(filtered_model_list)){
  model_posts[[i]] = as_draws_df(filtered_model_list[[i]]) %>% 
    mutate(model = names(filtered_model_list[i]),
           weight = filtered_model_list[[i]]$data2$model_weight)
}

post_average_parameters = bind_rows(model_posts) %>% 
  # filter(.draw <= 1000) %>% 
  select(.draw, starts_with("b_"), model, weight) %>% 
  mutate(across(starts_with("b_"), ~ replace_na(.x, 0))) %>%
  group_by(.draw) %>% 
  summarize(across(starts_with("b_"), ~ weighted.mean(.x, weight))) %>% 
  mutate(data = "fish + macros")

saveRDS(post_average_parameters, file = "posteriors/post_average_parameters.rds")

post_average_parameters_all = bind_rows(model_posts) %>% 
  # filter(.draw <= 1000) %>% 
  select(.draw, starts_with("b_"), starts_with("sd_"), starts_with("r_"), model, weight) %>% 
  mutate(across(starts_with("b_"), ~ replace_na(.x, 0))) %>%
  mutate(across(starts_with("sd_"), ~ replace_na(.x, 0))) %>% 
  mutate(across(starts_with("r_"), ~ replace_na(.x, 0))) %>%
  group_by(.draw) %>% 
  summarize(across(starts_with("b_"), ~ weighted.mean(.x, weight)),
            across(starts_with("sd_"), ~ weighted.mean(.x, weight)),
            across(starts_with("r_"), ~ weighted.mean(.x, weight))) %>% 
  mutate(data = "fish + macros")


saveRDS(post_average_parameters_all, file = "posteriors/post_average_parameters_all.rds")

# plot 
post_average_parameters %>%
  select(-b_Intercept) %>% 
  pivot_longer(starts_with("b_")) %>% 
  group_by(name) %>% 
  mutate(median_value = median(value), na.rm = T) %>% 
  ggplot(aes(y = reorder(name, median_value), x = value, color = median_value)) + 
  stat_pointinterval() +
  geom_vline(xintercept = 0)

