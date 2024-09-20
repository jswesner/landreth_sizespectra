brm_dummy = readRDS("models/brm_landreth_fishmacros_topo.rds")

library(rethinking)
# simulate data

x1 = rparetocounts(lambda = -1.8) # `lambda` is required wording from brms. in this case it means the lambda exponent of the ISD
x2 = rparetocounts(lambda = -1.5)
x3 = rparetocounts(lambda = -1.2)
x4 = rparetocounts(lambda = -0.9)
x5 = rparetocounts(lambda = -0.6)

isd_data = tibble(x1 = x1,
                  x2 = x2,
                  x3 = x3,
                  x4 = x4,
                  x5 = x5) |> 
  pivot_longer(cols = everything(), names_to = "group", values_to = "x") |> 
  group_by(group) |> 
  mutate(xmin = min(x),
         xmax = max(x)) |> 
  group_by(group, x) |> 
  add_count(name = "counts") %>% 
  mutate(cont_pred = parse_number(group))

brm_test = update(brm_test, newdata = isd_data)
brm_test_cont = update(brm_test, newdata = isd_data, formula = . ~ cont_pred)

brm_test_int = update(brm_test_int, newdata = isd_data, 
                  formula = . ~ 1)

newdata = isd_data %>% ungroup %>% 
  sample_n(size = 1000, weight = counts, replace = T) %>% 
  mutate(counts = 1)

brm_list = list(brm_test, brm_test_int, brm_test_cont)

for(i in 1:length(brm_list)){
  brm_list[[i]]$criteria = c(waic(log_lik(brm_list[[i]], newdata = newdata, ndraws = 500)))
}

brm_list[[1]]$criteria$estimates
brm_list[[2]]$criteria$estimates
brm_list[[3]]$criteria$estimates

test_posts = newdata %>% 
  mutate(counts = 1) %>% 
  distinct(group, xmin, xmax, counts) %>% 
  add_epred_draws(brm_test) %>% 
  ungroup %>% 
  select(.epred, .draw, group)

logprob = newdata %>% ungroup %>% 
  mutate(id = 1:nrow(.)) %>% 
  left_join(test_posts %>% filter(.draw <= 500)) %>% 
  rowwise() %>% 
  mutate(log_lik = dparetocounts(x, lambda = .epred, xmin = xmin, xmax = xmax))

n_samples = 500

lppd = logprob %>% group_by(x, id) %>% 
  reframe(lppd = log_sum_exp(log_lik) - n_samples)

pwaic = logprob %>% group_by(x, id) %>% 
  reframe(pwaic = var(log_lik))

-2*sum(lppd$lppd - sum(pwaic$pwaic))





# recapture results of McElreath pg 222 -----------------------------------

library(rethinking)
data(cars)

m <- quap(
  alist(
    dist ~ dnorm(mu, sigma),
    mu <- a + b*speed,
    a ~ dnorm(0, 100),
    b ~ dnorm(0, 10),
    sigma ~ dexp(1)
  ), data = cars )

set.seed(94)

post = extract.samples(m, n = 1000)


logprob_222_list = sapply(1:n_samples,
                 function(s){
                   mu = post$a[s] + post$b[s]*cars$speed
                   dnorm(cars$dist, mu, post$sigma[s], log = TRUE)
                 })


logprob_222 = as_tibble(logprob_222_list) %>% 
  mutate(dist = cars$dist) %>% 
  pivot_longer(cols = -dist) %>% 
  mutate(.draw = parse_number(name),
         method = "mcelreath_base",
         log_lik = value)

lppd = sapply(1:nrow(cars), function(i) log_sum_exp(logprob_222_list[i,]) - log(n_samples))
pwaic = sapply(1:nrow(cars), function(i) var(logprob_222_list[i,]))

-2*(sum(lppd) - sum(pwaic))


# Confirm that the tidy method gets the same answer -----------------------------

logprob = as_tibble(post) %>% 
  mutate(.draw = 1:nrow(.)) %>% 
  cross_join(cars %>% mutate(id = 1:nrow(.))) %>% 
  mutate(mu = a + b*speed,
         log_lik = dnorm(dist, mean = mu, sd = sigma, log = T),
         method = "mcelreath_tidy")


n_samples = max(logprob$.draw)

lppd_tidy = logprob %>% 
  group_by(dist, id) %>% 
  reframe(lppd = log_sum_exp(log_lik) - log(n_samples))

pwaic_tidy = logprob %>% 
  group_by(dist, id) %>% 
  reframe(pwaic = var(log_lik))

-2*(sum(lppd_tidy$lppd) - sum(pwaic_tidy$pwaic))
-2*(sum(lppd) - sum(pwaic))



# Get waic of isd model using tidy method ---------------------------------

get_waic_intonly = function(model, nsamples = 10, ndraws = 500){ 
  newdata = model$data %>% ungroup%>% 
    sample_n(size = nsamples, weight = counts, replace = T) %>% 
    mutate(counts = 1)
  
  post = as_draws_df(model) %>% filter(.draw <= ndraws)
  
  logprob = cross_join(post, newdata) %>% 
    rename(.epred = b_Intercept) %>% 
    rowwise() %>% 
    mutate(log_lik = dparetocounts(x, lambda = .epred, xmin = xmin, xmax = xmax),
           method = "mcelreath_tidy")
  
  n_samples = max(logprob$.draw)
  
  lppd = logprob %>% 
    group_by(x) %>% 
    reframe(lppd = log_sum_exp(log_lik) - log(n_samples))
  
  pwaic = logprob %>% 
    group_by(x) %>% 
    reframe(pwaic = var(log_lik))
  
  waic = -2*(sum(lppd$lppd) - sum(pwaic$pwaic)) 
  waic_se = sqrt(nrow(newdata)*var(-2*lppd$lppd - pwaic$pwaic))
  isd_waic = tibble(waic = waic, waic_se = waic_se)
  
  rm(logprob, n_samples, lppd, pwaic, newdata, waic, waic_se)
  
  print(isd_waic)
  
}

get_waic_isd = function(model, nsamples = 1000, ndraws = 500){
  newdata = model$data %>% ungroup%>% 
    sample_n(size = nsamples, weight = counts, replace = T) %>% 
    mutate(counts = 1)
  
  logprob = newdata %>% rownames_to_column(var = "id") %>% 
    add_epred_draws(model, ndraws = ndraws) %>% 
    rowwise() %>% 
    mutate(log_lik = dparetocounts(x, lambda = .epred, xmin = xmin, xmax = xmax),
           method = "mcelreath_tidy")
  
  n_samples = max(logprob$.draw)
  
  lppd = logprob %>% 
    group_by(x) %>% 
    reframe(lppd = log_sum_exp(log_lik) - log(n_samples))
  
  pwaic = logprob %>% 
    group_by(x) %>% 
    reframe(pwaic = var(log_lik))
  
  waic = -2*(sum(lppd$lppd) - sum(pwaic$pwaic)) 
  waic_se = sqrt(nrow(newdata)*var(-2*lppd$lppd - pwaic$pwaic))
  isd_waic = tibble(waic = waic, waic_se = waic_se)
  
  rm(logprob, n_samples, lppd, pwaic, newdata, waic, waic_se)
  
  print(isd_waic)
}

get_waic_intonly(brm_test_int, ndraws = 5)
get_waic_isd(brm_test, nsamples = 1000, ndraws = 5)
get_waic_isd(brm_test_cont, nsamples = 1000, ndraws = 500)


