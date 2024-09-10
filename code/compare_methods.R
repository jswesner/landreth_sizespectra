library(sizeSpectra)
library(tidyverse)
library(janitor)
library(tidybayes)

# load data
landreth_fishmacros_data = readRDS(file = "data/landreth_fishmacros_data.rds") 

# Get isdbayes results ----------------------------------------------------
#format post dots to fit sizespectra outputs
post_dots_sizespectra = readRDS(file = "posteriors/post_dots.rds") %>%
  filter(grepl("ntercept", model_name)) %>% 
  filter(taxon == "fish + macros") %>% 
  group_by(stream) %>% 
  median_qi(.epred) %>% 
  mutate(Year = 1980,
         Method = "isdbayes",
         b = .epred,
         confMin = .lower,
         confMax = .upper) %>% 
    select(Year, Method, b, confMin, confMax, stream)


# Get eight.methods results (Edwards 2017) --------------------------------
sizespectra_dat_temp = landreth_fishmacros_data %>%
  group_by(stream) %>%
  sample_n(500, weight = counts, replace = T) 

sizespectra_dat = sizespectra_dat_temp %>% 
  mutate(bodyMass = dw_g,
         Number = 1,
         Year = 1980) %>% 
  select(bodyMass, Year, stream, Number) %>% 
  group_by(stream) %>% 
  group_split()

eight_results = NULL
for(i in 1:length(sizespectra_dat)){
  eight_results[[i]] = eightMethods.count(sizespectra_dat[[i]]) %>% 
    as_tibble() %>% 
    mutate(stream = unique(sizespectra_dat[[i]]$stream))
}

saveRDS(eight_results, file = "data/eight_results.rds")

# Get mle.bin results (Edwards 2020 using Justin's code) --------------------------------------
dat_for_mlebin = sizespectra_dat_temp %>% 
  select(dw_g, stream) %>% 
  mutate(count = 1) %>% 
  group_by(stream) %>% 
  group_split() 

mle_bin_results = NULL

for(i in 1:length(dat_for_mlebin)){
  x_binned = binData(x = dat_for_mlebin[[i]]$dw_g,
                     binWidth = "2k")
  
  num_bins <- nrow(x_binned$binVals)
  
  # bin breaks are the minima plus the max of the final bin:
  bin_breaks <- c(dplyr::pull(x_binned$binVals, binMin),
                  dplyr::pull(x_binned$binVals, binMax)[num_bins])
  
  bin_counts <- dplyr::pull(x_binned$binVals, binCount)
  
  mle_temp_bin <-  calcLike(negLL.PLB.binned, # this is basically the mleBin part
                            p = -1.5,
                            w = bin_breaks, # calculated above
                            d = bin_counts, # calculated above
                            J = length(bin_counts),   # = num.bins
                            vecDiff = 1, # increase this if hit a bound, I think in aquasync I set it to ~10
                            suppress.warnings = TRUE) 
  
  mle_bin_results[[i]] = tibble(b = mle_temp_bin$MLE,
                                 confMin = mle_temp_bin$conf[1],
                                 confMax = mle_temp_bin$conf[2],
                                 Method = "mle.bin",
                                 stream = unique(dat_for_mlebin[[i]]$stream))
}

# get Landreth original log-log bin results ----------------------------------------------------
landreth_slopes = read_csv("data/landreth_slopes.csv") %>% 
  select(stream, macros_slope) %>% 
  rename(b = macros_slope) %>% 
  mutate(b = b - 1) %>% 
  mutate(Method = "Landreth tbl6 - 1")

# Combine and compare -----------------------------------------------------
all_predictors = read_csv("data/all_predictors.csv")

sizespectra_compare = bind_rows(eight_results) %>% 
  bind_rows(post_dots_sizespectra) %>% 
  bind_rows(bind_rows(mle_bin_results)) %>% 
  bind_rows(landreth_slopes) %>% 
  # filter(Method %in% c("isdbayes", "MLE", "Landreth original", "mle.bin")) %>%
  filter(Method != "LT") %>% 
  filter(Method != "LTplus1") %>% 
  filter(Method != "Llin") %>% 
  filter(Method != "LBmiz") 
  
compare_methods_plot= sizespectra_compare %>% 
  left_join(all_predictors) %>% 
  pivot_longer(cols = contains("PC")) %>%  
  mutate(Method = as.factor(Method),
         Method = fct_relevel(Method, "Landreth tbl6 - 1", "isdbayes", "MLE", "mle.bin")) %>% 
  ggplot(aes(x = value, y = b, ymin = confMin, ymax = confMax, group = Method, fill = Method, color = Method)) +
  geom_point(position = position_dodge(width = 0.5), size = 0.2, color = "grey50") +
  ggh4x::facet_grid2(name~Method) +
  geom_smooth(method = lm, linewidth = 0.1) + 
  brms::theme_default() +
  theme(strip.text = element_text(size = 7)) + 
  labs(y = "\u03bb", 
       x = "Predictor value (z-score)",
       caption = "Figure SX. Simple linear regression between lambda (i.e., the size spectra slope) and each PCA variable. Lambdas are produced by 7 methods,
       including four binning methods and 3 non-binning methods (isdbayes, MLE, mle.bin). The variation produced by different 
       methods leads to variation in the implied relationships of lambda with predictor variables.") +
  guides(fill = "none", 
         color = "none")

ggsave(compare_methods_plot, file = "plots/compare_methods_plot.jpg", width = 6.5, height = 9)

compare_methods_lambdas = sizespectra_compare %>% 
  select(stream, Method, b) %>% 
  pivot_wider(names_from = Method, values_from = b) %>% 
  pivot_longer(cols = c(LBbiom, LBNbiom, LCD, MLE, `Landreth tbl6 - 1`, mle.bin),
               names_to = "Method", values_to = "b") %>% 
  ggplot(aes(x = isdbayes, y = b, color = Method)) + 
  geom_point() + 
  facet_wrap(~Method) +
  geom_smooth(method = lm) +
  labs(y = "\u03bb",
       caption = "Figure SX. Correlations of CSS slopes produced by different methods compared to our Bayesian isdbayes method. All are positively correlated, but there is also
       substantial variation in the lambda produced by a given method. The best correlations are for the LCD, MLE, and mle.bin, 
       since these are all versions of a non-binning approach. ") +
  brms::theme_default() +
  guides(fill = "none", 
         color = "none")

ggsave(compare_methods_lambdas, file = "plots/compare_methods_lambdas.jpg")

# run models with landreth ------------------------------------------------

# get Landreth original log-log bin results ----------------------------------------------------
all_predictors = read_csv("data/all_predictors.csv")

landreth_slopes = read_csv("data/landreth_slopes.csv") %>% 
  select(stream, macros_slope) %>% 
  rename(b = macros_slope) %>% 
  mutate(b = b - 1) %>% 
  mutate(Method = "Landreth tbl6 - 1") %>% 
  left_join(all_predictors)

brm_landreth_herbivore = brm(b ~ fishmacros_herbivore_s, data = landreth_slopes)
brm_landreth_topo = update(brm_landreth_herbivore, formula = . ~ p_mda_elev_s + p_mda_slope_s + p_da_s,
                           newdata = landreth_slopes)
brm_landreth_invertivore = update(brm_landreth_herbivore, formula = . ~ fishmacros_invertivore_s,
                           newdata = landreth_slopes)
brm_landreth_elev_std = update(brm_landreth_herbivore, formula = . ~ p_elva_std_s ,
                                  newdata = landreth_slopes)
brm_landreth_landuse = update(brm_landreth_herbivore, formula = . ~ p_devel_s + p_ag_s ,
                               newdata = landreth_slopes)
brm_landreth_intercept = update(brm_landreth_herbivore, formula = . ~ 1 ,
                              newdata = landreth_slopes)
brm_landreth_watershedpca = update(brm_landreth_herbivore, formula = . ~ WatershedPC1_s + WatershedPC2_s + WatershedPC3_s + WatershedPC4_s + WatershedPC5_s ,
                              newdata = landreth_slopes)
brm_landreth_omnivore = update(brm_landreth_herbivore, formula = . ~ fishmacros_omnivore_s, 
                               newdata = landreth_slopes)
brm_landreth_slope_std = update(brm_landreth_herbivore, formula = . ~ p_slope_std_s + p_mda_slope_s, 
                               newdata = landreth_slopes)
brm_landreth_watershed1 = update(brm_landreth_herbivore, formula = . ~ WatershedPC1_s, 
                                newdata = landreth_slopes)
brm_landreth_biota = update(brm_landreth_herbivore, formula = . ~ CombinedBiotiaPC1_s, 
                                 newdata = landreth_slopes)
brm_landreth_biota2 = update(brm_landreth_herbivore, formula = . ~ CombinedBiotiaPC2_s, 
                            newdata = landreth_slopes)
brm_landreth_biota1234 = update(brm_landreth_herbivore, formula = . ~ CombinedBiotiaPC1_s +
                                  CombinedBiotiaPC2_s + CombinedBiotiaPC3_s + CombinedBiotiaPC4_s, 
                             newdata = landreth_slopes)


landreth_models = list(brm_landreth_herbivore, brm_landreth_topo, brm_landreth_invertivore,
                       brm_landreth_elev_std, brm_landreth_landuse, brm_landreth_intercept,
                       brm_landreth_watershedpca, brm_landreth_omnivore, brm_landreth_slope_std,
                       brm_landreth_watershed1, brm_landreth_biota, brm_landreth_biota2,
                       brm_landreth_biota1234)

saveRDS(landreth_models, file = "models/temporary/landreth_models.rds")

landreth_waic = NULL

for(i in 1:length(landreth_models)){
  temp_waics = waic(landreth_models[[i]])
  landreth_waic[[i]] = as_tibble(temp_waics$estimates) %>% 
    mutate(method = c("elpd_waic", "pwaic", "waic"),
           model = as.character(landreth_models[[i]]$formula)[1])
}

bind_rows(landreth_waic) %>% filter(method == "waic") %>% arrange(Estimate) %>% 
  ggplot(aes(x = reorder(model, Estimate), y = Estimate, ymin = Estimate - SE, ymax = Estimate + SE)) + 
  geom_pointrange() +
  coord_flip()

# check with other methods ----------------------------------------
eight_results_bind = bind_rows(readRDS(file = "data/eight_results.rds")) %>% 
  left_join(all_predictors)

####lbnbiom
lbnbiom = eight_results_bind %>% filter(Method == "LBNbiom")

lbnbiom_models = NULL

for(i in 1:length(landreth_models)){
  lbnbiom_models[[i]] = update(landreth_models[[i]], newdata = lbnbiom)
}

lbnbiom_waic = NULL

for(i in 1:length(lbnbiom_models)){
  temp_waics = waic(lbnbiom_models[[i]])
  lbnbiom_waic[[i]] = as_tibble(temp_waics$estimates) %>% 
    mutate(method = c("elpd_waic", "pwaic", "waic"),
           model = as.character(lbnbiom_models[[i]]$formula)[1])
}


####lbbiom
lbbiom = eight_results_bind %>% filter(Method == "LBbiom")

lbbiom_models = NULL

for(i in 1:length(landreth_models)){
  lbbiom_models[[i]] = update(landreth_models[[i]], newdata = lbbiom)
}

lbbiom_waic = NULL

for(i in 1:length(lbbiom_models)){
  temp_waics = waic(lbbiom_models[[i]])
  lbbiom_waic[[i]] = as_tibble(temp_waics$estimates) %>% 
    mutate(method = c("elpd_waic", "pwaic", "waic"),
           model = as.character(lbbiom_models[[i]]$formula)[1])
}



###MLE
mle = eight_results_bind %>% filter(Method == "MLE")

mle_models = NULL

for(i in 1:length(landreth_models)){
  mle_models[[i]] = update(landreth_models[[i]], newdata = mle)
}

mle_waic = NULL

for(i in 1:length(mle_models)){
  temp_waics = waic(mle_models[[i]])
  mle_waic[[i]] = as_tibble(temp_waics$estimates) %>% 
    mutate(method = c("elpd_waic", "pwaic", "waic"),
           model = as.character(mle_models[[i]]$formula)[1])
}



###isd

isd = post_dots_sizespectra %>% 
  left_join(all_predictors)

isd_models = NULL

for(i in 1:length(landreth_models)){
  isd_models[[i]] = update(landreth_models[[i]], newdata = isd)
}

isd_waic = NULL

for(i in 1:length(isd_models)){
  temp_waics = waic(isd_models[[i]])
  isd_waic[[i]] = as_tibble(temp_waics$estimates) %>% 
    mutate(method = c("elpd_waic", "pwaic", "waic"),
           model = as.character(isd_models[[i]]$formula)[1])
}


####combine
all_waics = bind_rows(bind_rows(landreth_waic) %>% mutate(b_method = "landreth"),
          bind_rows(isd_waic) %>% mutate(b_method = "isd_twostep"),
          bind_rows(mle_waic) %>% mutate(b_method = "mle"),
          bind_rows(lbbiom_waic) %>% mutate(b_method = "lbbiom"),
          bind_rows(lbnbiom_waic) %>% mutate(b_method = "lbnbiom")) %>% 
  filter(method == "waic") 
  


# show top two models for each method
fishmacros_waic_tocompare = readRDS(file = "posteriors/fishmacros_waic.rds") %>% 
  rename(Estimate = waic,
         SE = waic_se, 
         model = model_formula) %>% 
  select(Estimate, SE, model) %>% 
  mutate(b_method = "isdbayes") %>%
  mutate(model = str_replace(model, "dw_g \\| vreal\\(counts, xmin, xmax\\) ", "b "),
         model = str_remove(model, "\\+ \\(1 \\| stream\\) \\+ \\(1 \\| watershed\\)"),
         model = str_replace(model, " \\(1 \\| stream\\) \\+ \\(1 \\| watershed\\)", " 1"),
         model = str_trim(model),
         model = as.factor(model))


all_waics_df = bind_rows(all_waics, fishmacros_waic_tocompare)
saveRDS(all_waics_df, file = "models/temporary/all_waics_df.rds")  

all_waics_df %>% 
  group_by(b_method) %>% 
  arrange(b_method, Estimate) %>% 
  mutate(model_order = row_number()) %>% 
  filter(model_order <= 2) %>% 
  select(model_order, b_method, model) 


all_waics_df %>% 
  ggplot(aes(x = reorder(model, Estimate), 
             y = Estimate, ymin = Estimate - SE, ymax = Estimate + SE,
             color = b_method)) + 
  geom_pointrange() +
  coord_flip() +
  facet_wrap(~b_method, scales = "free_x") +
  guides(color = "none")

# check weighted regressions (need to finish)
brm_lbn_elev_std = update(brm_landreth_herbivore, formula = . ~ p_elva_std_s ,
                               newdata = lbnbiom)
brm_lbn_elev_std_mi = update(brm_landreth_herbivore, formula = b|mi(stdErr) ~ p_elva_std_s ,
                          newdata = lbnbiom)

plot(conditional_effects(brm_lbn_elev_std_mi), points = T)

