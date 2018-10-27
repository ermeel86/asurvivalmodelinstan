library(pec)
library(HSAUR)
library(tibble)
library(dplyr)
library(rstan)
library(purrr)
library(survival)
library(bayesplot)
rstan_options(auto_write = TRUE)
###############################################################################
data("mastectomy")
df <- as.tibble(mastectomy)
df <- df %>% mutate(metastized=as.double(metastized=="yes"))
###############################################################################
sm <- stan_model("~/Desktop/Stan/A_Survival_Model_in_Stan/exponential_survival_simple_ppc_pec.stan")
sm2 <- stan_model("~/Desktop/Stan/A_Survival_Model_in_Stan/weibull_survival_ppc_pec.stan")
sm3 <- stan_model("~/Desktop/Stan/A_Survival_Model_in_Stan/gamma_survival_ppc_pec.stan")
###############################################################################
N <- nrow(df)
X <- as.matrix(pull(df, metastized))
is_censored <- pull(df,event)==0
times <- pull(df,time)
msk_censored <- is_censored == 1
N_censored <- sum(msk_censored)
unique_times <- c(0,sort(unique(times)))
###############################################################################
stan_data <- list(N_uncensored=N-N_censored, 
                  N_censored=N_censored, 
                  X_censored=as.matrix(X[msk_censored,]),
                  X_uncensored=as.matrix(X[!msk_censored,]),
                  times_censored=times[msk_censored],
                  times_uncensored = times[!msk_censored],
                  NC=ncol(X),
                  N_times_eval_pec = length(unique_times),
                  times_eval_pec = unique_times
)
###############################################################################
# create ordered df so that it corresponds to Stan's internal ordering (to make pec work)
df_ordered <- tibble(time=c(stan_data$times_uncensored,stan_data$times_censored),
                     event=c(rep(1,N-N_censored), rep(0, N_censored)),
                     metastized=c(as.vector(X[!msk_censored,]), as.vector(X[msk_censored,]))
)
###############################################################################
fit <- sampling(sm, data=stan_data, seed=42, chains=4, cores=2, iter=4000)
fit2 <- sampling(sm2, data=stan_data, seed=42, chains=4, cores=2, iter=4000)
fit3 <- sampling(sm3, data=stan_data, seed=42, chains=4, cores=2, iter=4000)
###############################################################################
print(fit, pars=c("intercept", "betas[1]"))
print(fit2, pars=c("intercept", "betas[1]","alpha"))
print(fit3, pars=c("intercept", "betas[1]","alpha"))
###############################################################################
survs_var_names <- pull(mutate(as.tibble(expand.grid(x=1:length(unique_times) , y=1:nrow(X))),surv = sprintf("survs[%d,%d]",y,x)),"surv")
sprobs_mtx <- matrix(summary(fit)$summary[survs_var_names,"mean"], nrow=nrow(X), byrow=TRUE)
sprobs_mtx2 <- matrix(summary(fit2)$summary[survs_var_names,"mean"], nrow=nrow(X), byrow=TRUE)
sprobs_mtx3 <- matrix(summary(fit3)$summary[survs_var_names,"mean"], nrow=nrow(X), byrow=TRUE)
###############################################################################
rslt <- pec(object=list("expBayes"=sprobs_mtx,
                        "WeibullBayes"=sprobs_mtx2,
                        "GammaBayes"=sprobs_mtx3,
                        "expFreq"=coxph(Surv(time,event)~metastized,data=df_ordered,x=TRUE,y=TRUE)
                 ), 
            formula=Surv(time,event)~1, data=df_ordered, exact=TRUE, cens.model="marginal", 
            splitMethod="none",
            B=0,
            verbose=TRUE)
summary(rslt, times=seq(0, max(times), 5))
plot(rslt,xlim=c(0,200))
###############################################################################
surv_times_train <- times[!msk_censored]
post <- as.array(fit)
surv_times_rep1 <- as.matrix(map_dfr(1:dim(post)[2], ~as.tibble(post[,.,sprintf("times_uncensored_sampled[%d]", 1:stan_data$N_uncensored)])))
post <- as.array(fit2)
surv_times_rep2 <- as.matrix(map_dfr(1:dim(post)[2], ~as.tibble(post[,.,sprintf("times_uncensored_sampled[%d]", 1:stan_data$N_uncensored)])))
post <- as.array(fit3)
surv_times_rep3 <- as.matrix(map_dfr(1:dim(post)[2], ~as.tibble(post[,.,sprintf("times_uncensored_sampled[%d]", 1:stan_data$N_uncensored)])))
################################################################################
color_scheme_set("pink")
bayesplot::bayesplot_grid(
  plots = list(
    ppc_stat(surv_times_train, surv_times_rep1, binwidth = 1, stat = "mean"),
    ppc_stat(surv_times_train, surv_times_rep2, binwidth = 1, stat = "mean"),
    ppc_stat(surv_times_train, surv_times_rep3 , binwidth = 1, stat = "mean")
    ),
  titles = c("expBayes","WeibullBayes", "GammaBayes"),
  legends = FALSE,
  grid_args = list(ncol = 1)
)
################################################################################
bayesplot::bayesplot_grid(
  plots = list(
    ppc_stat(surv_times_train, surv_times_rep1, binwidth = 1, stat = "median"),
    ppc_stat(surv_times_train, surv_times_rep2, binwidth = 1, stat = "median"),
    ppc_stat(surv_times_train, surv_times_rep3 , binwidth = 1, stat = "median")
  ),
  titles = c("expBayes","WeibullBayes", "GammaBayes"),
  legends = FALSE,
  grid_args = list(ncol = 1)
)
################################################################################
bayesplot::bayesplot_grid(
  plots = list(
    ppc_stat(surv_times_train, surv_times_rep1, binwidth = 1, stat = "sd"),
    ppc_stat(surv_times_train, surv_times_rep2, binwidth = 1, stat = "sd"),
    ppc_stat(surv_times_train, surv_times_rep3 , binwidth = 1, stat = "sd")
  ),
  titles = c("expBayes","WeibullBayes", "GammaBayes"),
  legends = FALSE,
  grid_args = list(ncol = 1)
)
################################################################################
bayesplot::bayesplot_grid(
  plots = list(
    ppc_stat(surv_times_train, surv_times_rep1, binwidth = 1, stat = "min"),
    ppc_stat(surv_times_train, surv_times_rep2, binwidth = 1, stat = "min"),
    ppc_stat(surv_times_train, surv_times_rep3 , binwidth = 1, stat = "min")
  ),
  titles = c("expBayes","WeibullBayes", "GammaBayes"),
  legends = FALSE,
  grid_args = list(ncol = 1)
)
################################################################################
bayesplot::bayesplot_grid(
  plots = list(
    ppc_stat(surv_times_train, surv_times_rep1, binwidth = 1, stat = "max"),
    ppc_stat(surv_times_train, surv_times_rep2, binwidth = 1, stat = "max"),
    ppc_stat(surv_times_train, surv_times_rep3 , binwidth = 1, stat = "max")
  ),
  titles = c("expBayes","WeibullBayes", "GammaBayes"),
  legends = FALSE,
  grid_args = list(ncol = 1)
)
################################################################################