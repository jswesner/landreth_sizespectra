library(tidybayes)
library(tidyverse)
library(brms)

source("code/custom_functions.R")

# load data and models -----------------------------------------
# load data
landreth_fishmacros_data = readRDS(file = "data/landreth_fishmacros_data.rds") %>% mutate(taxon = "fish + macros")
landreth_macros_data = readRDS(file = "data/landreth_macros_data.rds") %>% mutate(taxon = "macros")
landreth_fish_data = readRDS(file = "data/landreth_fish_data.rds") %>% mutate(taxon = "fish")
landreth_data = bind_rows(landreth_fishmacros_data, landreth_macros_data, landreth_fish_data)

model_list_fishmacros = readRDS(file = "models/temporary/filtered_model_list.rds")
model_list_fish = readRDS(file = "models/fish_models/temporary/filtered_model_fish_list.rds")
model_list_macros = readRDS(file = "models/macros_models/temporary/filtered_model_macros_list.rds")

model_list = c(model_list_fish, model_list_macros, model_list_fishmacros)

# get stream lambdas ------------------------------------------------------
brm_dots_list = NULL

for(i in 1:length(model_list)){
  brm_dots_list[[i]] = model_list[[i]]$data %>% 
    select(-dw_g, -counts) %>% 
    mutate(model = as.character(model_list[[i]]$formula)[[1]],
           taxon = as.character(names(model_list[i]))) %>% 
    distinct() %>% 
    mutate(counts = 1)
}

epred_list = NULL
for(i in 1:length(model_list)){
  model_name = names(model_list[i])
    
  epred_list[[i]] = brm_dots_list[[i]] %>% add_epred_draws(model_list[[i]], re_formula = NULL) %>% 
    mutate(model_number = i,
           model_name = model_name) %>% 
    ungroup
}

post_dots = bind_rows(epred_list) %>% as_tibble() %>%
  mutate(intercept_s = case_when(grepl("ntercept", model_name) ~ 0)) %>% # dummy column to allow intercept only model 
  pivot_longer(cols = ends_with("_s")) %>% 
  filter(!is.na(value)) %>% 
  mutate(name = as.factor(name),
         name = fct_relevel(name, "intercept_s")) %>% 
  mutate(model = str_remove(model, "dw_g \\| vreal\\(counts, xmin, xmax\\) "),
         model = str_replace(model, "\\+ \\(1 \\| stream\\) \\+ \\(1 \\| watershed\\)", "..."),
         model = str_replace(model, " \\(1 \\| stream\\) \\+ \\(1 \\| watershed\\)", "Intercept"),
         model = as.factor(model), 
         model = fct_relevel(model, "~Intercept")) 

saveRDS(post_dots, file = "posteriors/post_dots.rds")

# get stream regression posteriors ------------------------
# 
# lines_list = NULL
# for(i in c(1, 2, 4, 5, 6)){
#   lines_init = conditional_effects(model_list[[i]])[1]
#   lines_list[[i]] = lines_init[[1]] %>% as_tibble() %>% 
#     mutate(model = as.character(model_list[[i]]$formula)[[1]]) %>% 
#     pivot_longer(cols = ends_with("_s")) 
# }

lines_list <- list()  

model_list_noint = model_list[!grepl("intercept", names(model_list))]

for (j in 1:length(model_list_noint)) {
  # Get the conditional effects for the current model
  lines_init <- conditional_effects(model_list_noint[[j]])
  
  # Iterate over each effect in lines_init
  for (k in 1:length(lines_init)) {  
    # Get the name of the current effect
    first_col_name <- names(lines_init)[k]
    
    # Add the current effect to lines_list with additional model information
    lines_list[[length(lines_list) + 1]] <- lines_init[[k]] %>%
      as_tibble() %>%
      mutate(model = as.character(model_list_noint[[j]]$formula)[[1]],
             name = first_col_name)
  }
}


post_lines = bind_rows(lines_list) %>% 
  add_row(name = "intercept_s",
          model = "dw_g | vreal(counts, xmin, xmax) ~ (1 | stream) + (1 | watershed)") %>%
  mutate(name = as.factor(name),
         name = fct_relevel(name, "intercept_s")) %>%
  mutate(model = str_remove(model, "dw_g \\| vreal\\(counts, xmin, xmax\\) "),
         model = str_replace(model, "\\+ \\(1 \\| stream\\) \\+ \\(1 \\| watershed\\)", "..."),
         model = str_replace(model, " \\(1 \\| stream\\) \\+ \\(1 \\| watershed\\)", "Intercept"),
         model = as.factor(model),
         model = fct_relevel(model, "~Intercept"))

saveRDS(post_lines, file = "posteriors/post_lines.rds")

