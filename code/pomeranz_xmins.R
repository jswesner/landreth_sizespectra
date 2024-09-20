library(tidybayes)
library(tidyverse)
library(brms)
library(ggthemes)
library(ggh4x)

landreth_fishmacros_data = readRDS(file = "data/landreth_fishmacros_data.rds") %>% mutate(taxon = "fish + macros")
peak_min_sizes = readRDS("file/data/peak_min_sizes.rds") # count and sumnorm xmins from Pomeranz in Sep 2024. These are
# derived from the mlebin procedure in Edwards SizeSpectra. They identify the binned values above which the body sizes are most 
# likely to follow a power law.

landreth_fishmacros_data = readRDS(file = "data/landreth_fishmacros_data.rds") 



# compare log bins of original and modified data --------------------------

original_data = landreth_fishmacros_data %>% mutate(data = "original")
modified_data = landreth_fishmacros_data %>% 
  left_join(peak_min_sizes) %>% 
  filter(dw_g >= sum_min) %>%  # remove biased small sizes
  group_by(stream) %>% 
  mutate(xmin = min(dw_g),
         data = "modified")

bind_rows(original_data, modified_data) %>% 
  group_by(stream, data) %>% 
  sample_n(20000, replace = T, weight = counts) %>% 
  mutate(log2_bins = cut(log2(dw_g), breaks = 10)) %>% 
  group_by(log2_bins, stream, data)  %>%
  tally() %>% 
  separate(log2_bins, into = c("min", "max"), sep = ",", remove = F) %>% 
  mutate(min = 2^parse_number(min),
         max = 2^parse_number(max)) %>% 
  ggplot(aes(x = min, y = n, color = data)) + 
  facet_wrap(~stream) +
  geom_point(aes(shape = data)) +
  geom_line() +
  scale_x_log10() +
  scale_color_colorblind() +
  NULL

# compare lambdas of original and modified data ---------------------------

# 1) intercept only models
original_and_modified_data = bind_rows(original_data, modified_data)

newlist = original_and_modified_data %>% group_by(stream, data) %>% 
  mutate(xmin = min(dw_g)) %>% group_split()

brm_dummy = readRDS("models/temporary/brm_list.rds")[[1]]

brm_update = list()

for(i in 1:length(newlist)){
  brm_update[[i]] = update(brm_dummy, newdata = newlist[[i]],
                           data2 = list(stream = unique(newlist[[i]]$stream),
                                        data = unique(newlist[[i]]$data)))
}

saveRDS(brm_update, file = "models/temporary/brm_update.rds")


# 2) get posteriors
brm_compare_posts_list = list()

for(i in 1:length(brm_update)){
  brm_compare_posts_list[[i]] = brm_update[[i]]$data %>% ungroup %>% 
    distinct(xmin, xmax) %>% 
    mutate(counts = 1) %>% 
    add_epred_draws(brm_update[[i]]) %>% 
    mutate(data = brm_update[[i]]$data2$data,
           stream = brm_update[[i]]$data2$stream)
}

brm_compare_posts = bind_rows(brm_compare_posts_list)

#3) plot
brm_compare_posts %>% 
  ggplot(aes(x = stream, y = .epred, color = data)) + 
  stat_pointinterval()

brm_compare_posts %>% 
  group_by(stream, data) %>% 
  reframe(lambda = median(.epred)) %>% 
  pivot_wider(names_from = data, values_from = lambda) %>% 
  ggplot(aes(x = modified, y = original)) + 
  geom_point() +
  geom_abline() +
  ylim(-3, -1) + 
  xlim(-3, -1)

#4) check

dat_resample_list <- lapply(
  brm_update, function(x) {
  x$data %>%
    sample_n(1000, weight = counts, replace = TRUE) %>%
    mutate(stream = x$data2$stream,
           data = x$data2$data)
})

post_resample_list = list()

for(i in 1:length(dat_resample_list)){
  post_resample_list[[i]] = dat_resample_list[[i]] %>% 
    add_predicted_draws(brm_update[[i]], ndraws = 300)
}

bind_rows(post_resample_list) %>% 
  filter(.draw <= 10) %>% 
  ggplot(aes(x = .prediction)) + 
  geom_density(aes(group = .draw), alpha = 0.4) + 
  facet_grid2(stream ~ data, scales = "free_y") +
  scale_x_log10() +
  geom_density(data = bind_rows(dat_resample_list),
               aes(x = dw_g), color = "dodgerblue")

bind_rows(post_resample_list) %>% 
  group_by(stream, data, .draw) %>% 
  reframe(median = median(.prediction),
          gm = exp(mean(log(.prediction)))) %>% 
  ggplot(aes(x = gm)) + 
  facet_grid2(stream ~ data, scales = "free_y") +
  geom_histogram(bins = 30) +
  geom_vline(data = bind_rows(dat_resample_list) %>% group_by(stream, data) %>% reframe(median = median(dw_g),
                                                                                        gm = exp(mean(log(dw_g)))),
                 aes(xintercept = gm), color = "red3")
  
