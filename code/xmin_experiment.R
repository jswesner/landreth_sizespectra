# Load the packages
library(tidyverse)
library(brms)
library(isdbayes)

# analyze single site isd's with a site and global xmin, then compare
# result: lambda estimates do not depend on xmin choice

# load data
landreth_fishmacros_data = readRDS(file = "data/landreth_fishmacros_data.rds") 
landreth_fishmacros_data_globalxmin = landreth_fishmacros_data %>% 
  ungroup %>% 
  mutate(xmin = min(dw_g),
         xmin_method = "global")

dat = bind_rows(landreth_fishmacros_data %>% mutate(xmin_method = "persite"),
                landreth_fishmacros_data_globalxmin)

stream_list = dat %>% group_by(stream, xmin_method) %>% 
  group_split()

# fit (or load) dummy model
# brm_singles = brm(dw_g|vreal(counts, xmin, xmax) ~ 1,
#                   data = stream_list[[1]],
#                   stanvars = stanvars,
#                   family = paretocounts(),
#                   # prior = c(prior(normal(-1.2, 0.2), class = "Intercept"),
#                             # prior(exponential(5), class = "sd"),
#                             # prior(normal(0, 1), class = "b")),
#                   iter = 100, chains = 1,
#                   file = "models/temporary/brm_singles.rds")

brm_singles = readRDS("models/temporary/brm_singles.rds")

# update dummy model (avoids the need to recompile each time)
brm_singles_list = NULL

for(i in 1:length(stream_list)){
  brm_singles_list[[i]] = update(brm_singles, newdata = stream_list[[i]], 
                                 data2 = list(stream = unique(stream_list[[i]]$stream),
                                              xmin_method = unique(stream_list[[i]]$xmin_method)),
                                 iter = 500)
}

saveRDS(brm_singles_list, file = "models/temporary/brm_singles_list.rds")
brm_singles_list = readRDS(file = "models/temporary/brm_singles_list.rds")

# get posteriors
brm_singles_posts = NULL

for(i in 1:length(brm_singles_list)){
  brm_singles_posts[[i]] = as_draws_df(brm_singles_list[[i]]) %>% 
    mutate(stream = as.character(brm_singles_list[[i]]$data2[1]),
           xmin_method = as.character(brm_singles_list[[i]]$data2[2]),
           xmin_value = "per site")
}

brm_singles_posts_preds = bind_rows(brm_singles_posts)

brm_singles_posts_preds %>% 
  group_by(stream) %>% 
  mutate(median = median(b_Intercept)) %>% 
  ggplot(aes(x = reorder(stream, median), y = b_Intercept, color = xmin_method)) + 
  stat_pointinterval(position = position_dodge(width = 0.2))


