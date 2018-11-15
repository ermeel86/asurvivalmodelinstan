library(rstan);rstan_options(auto_write = TRUE)

sm <- stan_model("~/Desktop/Stan/A_Survival_Model_in_Stan/exponential_survival_cens_ppc.stan")

library(HSAUR)
library(tibble)
library(dplyr)
library(purrr)
library(bayesplot)
data("mastectomy")
df <- as.tibble(mastectomy)
df <- df %>% mutate(metastized=as.double(metastized=="yes"))

N <- nrow(df)
X <- as.matrix(pull(df, metastized))
is_censored <- pull(df,event)==0
times <- pull(df,time)
msk_censored <- is_censored == 1
N_censored <- sum(msk_censored)

stan_data <- list(N_uncensored=N-N_censored, 
                  N_censored=N_censored, 
                  X_censored=as.matrix(X[msk_censored,]),
                  X_uncensored=as.matrix(X[!msk_censored,]),
                  times_censored=times[msk_censored],
                  times_uncensored = times[!msk_censored],
                  NC=ncol(X)
)

fit <- sampling(sm, data=stan_data, seed=42, chains=4, cores=2, iter=4000)
post <- as.array(fit)

bayesplot::mcmc_intervals(post, pars=c("betas_c[1]","betas[1]"))

bayesplot::mcmc_intervals(post, pars=c("intercept_c","intercept"))

bayesplot::mcmc_intervals(post, pars=c("intercept_c","intercept"), transformations = function(x) exp(-x))



surv_times_rep <- as.matrix(map_dfr(1:dim(post)[2], ~as.tibble(post[,.,sprintf("times_samples[%d]", 1:N)])))
events_rep <- as.matrix(map_dfr(1:dim(post)[2], ~as.tibble(post[,.,sprintf("events[%d]", 1:N)])))


color_scheme_set("brightblue")
ppc_stat(times, surv_times_rep, binwidth = 1, stat = "mean")
ppc_stat(times, surv_times_rep, binwidth = 1, stat = "sd")
ppc_stat(times, surv_times_rep, binwidth = 1, stat = "max")
ppc_stat(times, surv_times_rep, binwidth = 1, stat = "min")
events <-  as.double(pull(df,event))
color_scheme_set("red")
ppc_stat(events, events_rep, binwidth = 1, stat = "mean")
ppc_stat(events, events_rep, binwidth = 1, stat = "sd")

