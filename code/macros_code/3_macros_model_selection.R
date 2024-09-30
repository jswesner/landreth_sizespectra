library(tidybayes)
library(tidyverse)
library(brms)
library(isdbayes)

source("code/custom_functions.R")

# load data and models -----------------------------------------
landreth_macros_data = readRDS(file = "data/landreth_macros_data.rds") 

model_macros_files = list.files(
  path = "models/macros_models/", 
  pattern = "\\.rds$", 
  full.names = TRUE
)

# Exclude "model_macros_list.rds". Otherwise the code below gives an error
model_macros_files = model_macros_files[basename(model_macros_files) != "model_macros_list.rds"]

model_macros_list = NULL

# Loop through each file and read the RDS into the list
for (file in model_macros_files) {
  model_macros_list[[file]] <- readRDS(file)
}



# get waic for intercept only

for(i in 1:length(model_macros_list)){
  model_macros_waic = get_waic_isd(model_macros_list[[i]], nsamples = 1000, ndraws = 500, resp = "dw_g") %>% 
    mutate(model_macros_name = names(model_macros_list[i]),
           model_macros_formula = as.character(model_macros_list[[i]]$formula)[1])
  
  model_macros_list[[i]]$criteria = model_macros_waic
}

saveRDS(model_macros_list, file = "posteriors/model_macros_list.rds")


# posterior averaging -----------------------------------------------------
model_macros_list = readRDS(file = "posteriors/model_macros_list.rds")

# filter to only subset of models
# filtered_model_macros_names <- names(model_macros_list)[!grepl("PC", names(model_macros_list), ignore.case = TRUE)]
filtered_model_macros_names <- names(model_macros_list)[grepl("agtopo|trophic|intercept", names(model_macros_list), ignore.case = TRUE)]
filtered_model_macros_list <- model_macros_list[filtered_model_macros_names]
saveRDS(filtered_model_macros_list, file = "models/macros_models/temporary/filtered_model_macros_list.rds")

# get model weights (adpated from https://github.com/paul-buerkner/brms/blob/master/R/model_weights.R)
waic_values = numeric(length(filtered_model_macros_list)) 

for(i in 1:length(filtered_model_macros_list)) waic_values[i] = filtered_model_macros_list[[i]]$criteria$waic

se_values = numeric(length(filtered_model_macros_list))
for(i in 1:length(filtered_model_macros_list)) se_values[i] = filtered_model_macros_list[[i]]$criteria$waic_se
min_waic <- min(waic_values)
adjusted_diffs <- (waic_values - min_waic) / se_values
weights <- exp(-adjusted_diffs / 2)
model_macros_weights <- weights / sum(weights)  # Normalize to sum to 1

waic_weights = tibble(waic = waic_values, model_macros_weights = model_macros_weights, model = names(filtered_model_macros_list))

# add model weights to data2 in models
for(i in 1:length(filtered_model_macros_list)) filtered_model_macros_list[[i]]$data2 = list(model_macros_weight = model_macros_weights[i])

model_macros_posts = list()

for(i in 1:length(filtered_model_macros_list)){
  model_macros_posts[[i]] = as_draws_df(filtered_model_macros_list[[i]]) %>% 
    mutate(model = names(filtered_model_macros_list[i]),
           weight = filtered_model_macros_list[[i]]$data2$model_macros_weight)
}

post_average_macros_parameters = bind_rows(model_macros_posts) %>% 
  # filter(.draw <= 1000) %>% 
  select(.draw, starts_with("b_"), model, weight) %>% 
  mutate(across(starts_with("b_"), ~ replace_na(.x, 0))) %>%
  group_by(.draw) %>% 
  summarize(across(starts_with("b_"), ~ weighted.mean(.x, weight))) %>% 
  mutate(data = "fish + macros")

saveRDS(post_average_macros_parameters, file = "posteriors/post_average_macros_parameters.rds")


# plot 
post_average_macros_parameters %>%
  select(-b_Intercept) %>% 
  pivot_longer(starts_with("b_")) %>% 
  group_by(name) %>% 
  mutate(median_value = median(value), na.rm = T) %>% 
  ggplot(aes(y = reorder(name, median_value), x = value, color = median_value)) + 
  stat_pointinterval() +
  geom_vline(xintercept = 0)

