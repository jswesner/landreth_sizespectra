library(tidybayes)
library(tidyverse)
library(brms)

source("code/custom_functions.R")

# load data and models -----------------------------------------
landreth_fishmacros_data = readRDS(file = "data/landreth_fishmacros_data.rds")
model_list = readRDS(file = "models/temporary/filtered_model_list.rds")

# posterior predictive checks ---------------------------------------------------------

# compare geometric mean
post_gm_list = NULL

for(i in 1:length(model_list)){
  post_gm_list[[i]] = isd_ppcheck(model_list[[i]], re_formula = NULL) %>% 
    mutate(model = as.character(model_list[[i]]$formula)[[1]])
}

post_gm = bind_rows(post_gm_list)

post_gm %>%
  ggplot(aes(x = gm)) +
  geom_histogram(data = . %>% 
                   filter(.draw > 0),
                 bins = 50) +
  scale_x_log10() +
  geom_vline(data = . %>% filter(.draw == 0),
             aes(xintercept = gm), color = 'red') +
  facet_grid(stream_watershed ~ model)

post_gm %>% 
  filter(.draw > 0) %>% 
  ggplot(aes(y = reorder(stream_watershed, median_gm), x = gm)) +
  stat_halfeye() +
  scale_x_log10() +
  geom_point(data = post_gm %>% filter(.draw == 0), color = 'red') +
  facet_wrap(~model)


# bayes_p 
post_gm %>% filter(.draw > 0) %>% 
  left_join(post_gm %>% filter(.draw == 0) %>% 
              rename(y_gm = gm) %>% 
              distinct(stream_watershed, y_gm, model)) %>%
  mutate(diff = gm - y_gm) %>%
  group_by(stream_watershed, model) %>% 
  reframe(bayes_p = sum(diff < 0)/max(.draw)) %>% 
  pivot_wider(names_from = model, values_from = bayes_p)

post_gm %>% filter(.draw > 0) %>% 
  left_join(post_gm %>% filter(.draw == 0) %>% 
              rename(y_gm = gm) %>% 
              distinct(stream_watershed, y_gm, model)) %>%
  mutate(diff = gm - y_gm) %>%
  group_by(stream_watershed, model) %>% 
  reframe(bayes_p = sum(diff < 0)/max(.draw)) %>% 
  group_by(model) %>% 
  median_qi(bayes_p)


# ppcheck by hand ---------------------------------------------------------
# load posterior averaged parameter values. fishmacros model only
post_average_parameters_all = readRDS(file = "posteriors/post_average_parameters_all.rds")
predictors_to_keep = tibble(predictors = names(post_average_parameters_all)) %>% filter(grepl("b_", predictors)) %>% 
  mutate(predictors = str_remove(predictors, "b_")) %>% 
  filter(predictors != "Intercept") %>% 
  pull()

landreth_fishmacros_predictors = readRDS(file = "data/landreth_fishmacros_data.rds") %>% 
  select(all_of(predictors_to_keep), stream, watershed) %>% 
  distinct() %>% 
  pivot_longer(cols = ends_with("_s"), values_to = "x_value", names_to = "predictors")

stream_offsets = post_average_parameters_all %>% 
  select(starts_with("r_stream"), .draw) %>% 
  pivot_longer(cols = starts_with("r_stream"), names_to = "stream", values_to = "stream_offset") %>% 
  mutate(stream = str_remove(stream, "r_stream\\["),
         stream = str_remove(stream, ",Intercept\\]"),
         stream = str_replace(stream, "\\.", " "))

watershed_offsets = post_average_parameters_all %>% 
  select(starts_with("r_watershed"), .draw) %>% 
  pivot_longer(cols = starts_with("r_watershed"), names_to = "watershed", values_to = "watershed_offset") %>% 
  mutate(watershed = str_remove(watershed, "r_watershed\\["),
         watershed = str_remove(watershed, ",Intercept\\]")) %>% 
  left_join(landreth_fishmacros_data %>% ungroup %>% distinct(watershed, stream), relationship = "many-to-many")


stream_watershed_offsets = watershed_offsets %>% left_join(stream_offsets)

post_average_lambdas = post_average_parameters_all %>% 
  filter(.draw <= 1000) %>% 
  select(starts_with("b_"), .draw) %>%
  rename(Intercept = b_Intercept) %>% 
  pivot_longer(cols = starts_with("b_"), names_to = "predictors", values_to = "parameter_value") %>% 
  mutate(predictors = str_remove(predictors, "b_")) %>% 
  left_join(landreth_fishmacros_predictors) %>% 
  left_join(stream_watershed_offsets) %>% 
  mutate(parameter_vector = parameter_value*x_value) %>% 
  select(-predictors, -parameter_value) %>% 
  group_by(Intercept, .draw, stream, watershed, watershed_offset, stream_offset) %>% 
  reframe(parameter_sum = sum(parameter_vector)) %>% 
  arrange(.draw) %>% 
  mutate(lambda = Intercept + parameter_sum + watershed_offset + stream_offset)


# pp_check global --------------------------------------------------------

sim_data = model_list[[2]]$data %>%
  select(stream, watershed, xmin, xmax, dw_g, counts) %>% 
  sample_n(size = 1000, replace = T, weight = counts) %>% 
  mutate(xmin = min(dw_g))

post_average_preds = post_average_lambdas %>% 
  left_join(sim_data, relationship = "many-to-many") %>% 
  rowwise() %>% 
  filter(.draw <= 100) %>% 
  mutate(dw_g = rparetocounts(1, lambda = lambda, xmin = xmin, xmax = xmax),
         data = "ypred")

# combine predicted data sets with the original resampled data set
post_average_all = post_average_preds %>% 
  select(.draw, stream, watershed, dw_g, data) %>% 
  bind_rows(sim_data %>% select(stream, watershed, dw_g) %>%
              mutate(data = "yraw",
                     .draw = 0))

post_average_all %>% 
  filter(.draw <= 10) %>% 
  ggplot(aes(x = dw_g)) + 
  geom_histogram(aes(group = .draw, color = data)) +
  scale_x_log10() +
  facet_wrap(~.draw)

post_summaries = post_average_all %>% 
  group_by(.draw) %>% 
  reframe(gm = exp(mean(log(dw_g))))
  
post_summaries %>% 
  filter(.draw > 0) %>% 
  ggplot(aes(x = gm)) +
  geom_histogram() +
  geom_vline(data = post_summaries %>% filter(.draw == 0),
             aes(xintercept = gm))

post_summaries %>% 
  filter(.draw > 0) %>% 
  mutate(yraw = post_summaries %>% filter(.draw == 0) %>% select(gm) %>% pull) %>% 
  mutate(diff = gm - yraw) %>% 
  reframe(bayes_p_global_gm = sum(diff > 0)/max(.draw))

# pp check by stream ------------------------------------------------------
sim_data_streams = model_list[[1]]$data %>%
  select(stream, watershed, xmin, xmax, dw_g, counts) %>% 
  group_by(stream) %>% 
  sample_n(size = 1000, replace = T, weight = counts) %>% 
  mutate(xmin = min(dw_g),
         counts = 1)

post_average_preds_stream = post_average_lambdas %>% 
  left_join(sim_data, relationship = "many-to-many") %>% 
  rowwise() %>% 
  filter(.draw <= 100) %>% 
  mutate(dw_g = rparetocounts(1, lambda = lambda, xmin = xmin, xmax = xmax),
         data = "ypred")

post_average_stream = post_average_preds_stream %>% 
  select(.draw, stream, watershed, dw_g, data) %>% 
  bind_rows(sim_data %>% select(stream, watershed, dw_g) %>%
              mutate(data = "yraw",
                     .draw = 0))

post_average_stream %>% 
  filter(.draw <= 10) %>% 
  ggplot(aes(x = .draw, y = dw_g, color = data)) + 
  geom_jitter(alpha = 1, width = 0.1, size = 0.2) + 
  scale_y_log10() +
  scale_color_brewer() +
  geom_boxplot(aes(group = .draw), width = 1, outlier.shape = NA) + 
  facet_wrap(~stream, ncol = 3, scales = "free_y") + 
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10),
                     labels = c("yraw", "2", "4", "6", "8", "10"))

# gm
post_gm_stream = post_average_stream %>% 
  group_by(stream, .draw) %>% 
  reframe(gm = exp(mean(log(dw_g))),
          # median = median(dw_g)
          )

post_gm_stream %>% 
  filter(.draw > 0) %>% 
  left_join(post_gm_stream %>% filter(.draw == 0) %>% 
              select(-.draw) %>% rename(yraw = gm)) %>% 
  mutate(diff = gm - yraw) %>% 
  group_by(stream) %>% 
  reframe(bayes_p = sum(diff>0)/max(.draw))


# # use mle_bin xmin --------------------------------------------------------
# mlebin_xmin = sim_data %>% 
#   group_by(stream) %>% 
#   mutate(log2_bin = cut(dw_g, breaks = 2^seq(floor(log2(min(dw_g))), ceiling(log2(max(dw_g))), by = 1))) %>%
#   group_by(log2_bin, stream) %>%
#   summarise(
#     count = n(),
#     bin_min = min(dw_g),
#     bin_max = max(dw_g)
#   ) %>%
#   arrange(desc(count)) %>% 
#   group_by(stream) %>% 
#   filter(count == max(count))
# 
# new_sims = 
#   sim_data %>% 
#   left_join(mlebin_xmin %>% select(stream, bin_min)) %>% 
#   filter(dw_g >= bin_min)
# 
# new_xmin_data = landreth_fishmacros_data %>% 
#   left_join(mlebin_xmin %>% select(stream, bin_min)) %>% 
#   filter(dw_g >= bin_min) %>% 
#   mutate(xmin = min(dw_g))
# 
# new_xmin_data_list = new_xmin_data %>% group_by(stream) %>% group_split()
# 
# brm_test = brm(dw_g | vreal(counts, xmin, xmax) ~ 1, 
#                data = new_xmin_data_list[[1]],
#                stanvars = stanvars,
#                family = paretocounts(),
#                iter = 1000, chains = 1)
# 
# brm_list = list() 
# for(i in 1:length(new_xmin_data_list)){
#   brm_list[[i]] = update(brm_test, newdata = new_xmin_data_list[[i]],
#                     data2 = list(stream = new_xmin_data_list[[i]]$stream,
#                                  xmin = "log2bin_max"))
# }
# 
# saveRDS(brm_list, file = "models/temporary/brm_list.rds")
# 
# brm_post_list = list()
# for(i in 1:length(brm_list)){
#   brm_post_list[[i]] = as_draws_df(brm_list[[i]]) %>% mutate(stream = unique(brm_list[[i]]$data2$stream)) 
# }
# 
# bind_rows(brm_post_list) %>% 
#   group_by(stream) %>% 
#   mutate(median = median(b_Intercept)) %>% 
#   ggplot(aes(x = reorder(stream, -median), y = b_Intercept)) +
#   stat_pointinterval() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 0))
# 
# # compare new xmin to original xmin lambdas
# post_dots_stream = post_average_lambdas %>% 
#   group_by(stream) %>% 
#   reframe(median = median(lambda),
#           lower = quantile(lambda, 0.025),
#           upper = quantile(lambda, 0.975)) %>% 
#   mutate(data = "original")
# 
# new_post_dots_stream = bind_rows(brm_post_list) %>% 
#   group_by(stream) %>% 
#   reframe(median = median(b_Intercept),
#          lower = quantile(b_Intercept, 0.025),
#          upper = quantile(b_Intercept, 0.975)) %>% 
#   mutate(data = "new")
# 
# 
# post_dots_all = bind_rows(post_dots_stream,
#           new_post_dots_stream) %>% 
#   left_join(landreth_fishmacros_data %>% ungroup %>% 
#               select(-dw_g, -counts,-xmin, -xmax, -n, -sample_area_m2) %>% 
#               distinct())
#   
# post_dots_all %>% 
#   ggplot(aes(x = reorder(stream, -median), y = median, ymin = lower, ymax = upper)) +
#   geom_pointrange(aes(color = data))
# 
# 
# post_dots_all %>% 
#   ggplot(aes(x = fishmacros_omnivore_s, y = median, ymin = lower, ymax = upper)) +
#   geom_pointrange(aes(color = data)) +
#   facet_wrap(~data) +
#   geom_smooth(method = lm)
# 
# 
# 
# # get new postpredicts
# new_xmin_postpredict = list()
# sim_newxmin = list()
# 
# for(i in 1:length(brm_list)){
# sim_newxmin[[i]] = brm_list[[i]]$data %>%
#   sample_n(size = 1000, replace = T, weight = counts) %>% 
#   mutate(xmin = min(dw_g),
#          model = i)
# 
# new_xmin_postpredict[[i]] = sim_newxmin[[i]] %>%
#   add_predicted_draws(brm_list[[i]]) %>% 
#   mutate(model = i)
# }
# 
# bind_rows(new_xmin_postpredict) %>% 
#   filter(.draw == as.integer(runif(1, 1, 500))) %>% 
#   ggplot(aes(x = .prediction)) + 
#   stat_histinterval() +
#   scale_x_log10() +
#   stat_histinterval(data = bind_rows(sim_newxmin), aes(x = dw_g), color = "dodgerblue",
#                fill = "dodgerblue",
#                alpha = 0.4) +
#   facet_wrap(~model)
# 
