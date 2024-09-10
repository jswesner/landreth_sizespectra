# from McElreath 2020 (page 222)

# isd_waic <- function(model, ndraws = 500) {
#   waic(log_lik(model, ndraws = ndraws))
#   }

# isd_waic <- function(model, ndraws = 500) {
#   log_lik_values <- log_lik(model, ndraws = ndraws)
#   
#   mean_log_lik <- apply(log_lik_values, 2, mean)
#   var_log_lik <- apply(log_lik_values, 2, var)
#   
#   lppd <- sum(mean_log_lik)
#   p_waic <- sum(var_log_lik)
#   waic <- -2 * (lppd - p_waic)
#   
#   n_cases <- nrow(model$data)
#   waic_vec <- -2 * (mean_log_lik - var_log_lik)
#   waic_se <- sqrt(n_cases * var(waic_vec))
#   
#   model_name <- deparse(substitute(model))
#   print(tibble(model = model_name, waic = waic, waic_se = waic_se))
#   
#   rm(lppd, p_waic, waic, n_cases, waic_vec, waic_se, log_lik_values, mean_log_lik, var_log_lik)
# }

# a function from rethinking
log_sum_exp <- function( x ) {
  xmax <- max(x)
  xsum <- sum( exp( x - xmax ) )
  xmax + log(xsum)
}

get_waic_intonly = function(model, resp = "x", nsamples = 1000, ndraws = 500){ 
  newdata = model$data %>% ungroup%>% 
    sample_n(size = nsamples, weight = counts, replace = T) %>% 
    mutate(counts = 1)
  
  post = as_draws_df(model) %>% filter(.draw <= ndraws)
  
  logprob = cross_join(post, newdata) %>% 
    rename(.epred = b_Intercept) %>% 
    rowwise() %>% 
    mutate(log_lik = dparetocounts(!!sym(resp), lambda = .epred, xmin = xmin, xmax = xmax),
           method = "mcelreath_tidy")
  
  n_samples = max(logprob$.draw)
  
  lppd = logprob %>% 
    group_by(!!sym(resp)) %>% 
    reframe(lppd = log_sum_exp(log_lik) - log(n_samples))
  
  pwaic = logprob %>% 
    group_by(!!sym(resp)) %>% 
    reframe(pwaic = var(log_lik))
  
  waic = -2*(sum(lppd$lppd) - sum(pwaic$pwaic)) 
  waic_se = sqrt(nrow(newdata)*var(-2*lppd$lppd - pwaic$pwaic))
  isd_waic = tibble(waic = waic, waic_se = waic_se)
  
  rm(logprob, n_samples, lppd, pwaic, newdata, waic, waic_se)
  
  return(isd_waic)
  
}

get_waic_isd = function(model, resp = "x", nsamples = 1000, ndraws = 500, re_formula = NULL){
  newdata = model$data %>% ungroup%>% 
    sample_n(size = nsamples, weight = counts, replace = T) %>% 
    mutate(counts = 1)
  
  logprob = newdata %>% rownames_to_column(var = "id") %>% 
    add_epred_draws(model, ndraws = ndraws, re_formula = re_formula) %>% 
    rowwise() %>% 
    mutate(log_lik = dparetocounts(x = !!sym(resp), lambda = .epred, xmin = xmin, xmax = xmax),
           method = "mcelreath_tidy")
  
  n_samples = max(logprob$.draw)
  
  lppd = logprob %>% 
    group_by(!!sym(resp)) %>% 
    reframe(lppd = log_sum_exp(log_lik) - log(n_samples))
  
  pwaic = logprob %>% 
    group_by(!!sym(resp)) %>% 
    reframe(pwaic = var(log_lik))
  
  waic = -2*(sum(lppd$lppd) - sum(pwaic$pwaic)) 
  waic_se = sqrt(nrow(newdata)*var(-2*lppd$lppd - pwaic$pwaic))
  isd_waic = tibble(waic = waic, waic_se = waic_se)
  
  rm(logprob, n_samples, lppd, pwaic, newdata, waic, waic_se)
  
  return(isd_waic)
}


isd_ppcheck = function(model, n_pred = 500, re_formula = NULL){
  model_name <- deparse(substitute(model))

  sample_post_ypred = model$data %>% as_tibble() %>% 
    select(-dw_g, -counts) %>% 
    distinct() %>% 
    mutate(counts = 1) %>% 
    add_predicted_draws(model, ndraws = n_pred, re_formula = re_formula) %>% 
    ungroup %>% 
    distinct(stream, watershed, .draw, .prediction) %>% 
    mutate(source = "ypred",
           model = model_name)
  
  resample_data = model$data %>%
    group_by(stream, watershed) %>% 
    sample_n(size = n_pred, replace = T, weight = counts) %>% 
    arrange(-dw_g) %>% 
    mutate(n_yx = 1:max(row_number()),
           n = max(row_number()),
           ismax = dw_g - xmax) %>% 
    rename(.prediction = dw_g) %>% 
    mutate(.draw = 0) %>% 
    select(stream, watershed, .draw, .prediction) %>% 
    mutate(source = "y",
           model = model_name)
  
  bind_rows(resample_data, sample_post_ypred) %>% 
    ungroup %>%
    group_by(.draw, stream, watershed) %>% 
    reframe(gm = exp(mean(log(.prediction))),
            median = median(.prediction)) %>% 
    mutate(stream_watershed = paste0(stream, "_", watershed)) %>% 
    group_by(stream_watershed) %>% 
    mutate(median_gm = median(gm),
           model = model_name)
  
}
