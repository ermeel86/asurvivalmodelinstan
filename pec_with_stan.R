library(pec)
library(HSAUR)
library(tibble)
library(dplyr)
library(rstan)
library(purrr)
library(survival)
library(bayesplot)
library(splines2)
library(cowplot)
rstan_options(auto_write = TRUE)
source("~/Documents/Talks/Bayesian_Survival_Models_in_Stan_Meetup_20181105/src/utils.R")
###############################################################################
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# To use for fills, add
scale_fill_manual(values=cbbPalette)
# To use for line and point colors, add
scale_colour_manual(values=cbbPalette)
###############################################################################
data("mastectomy")
df <- as.tibble(mastectomy)
df <- df %>% mutate(metastized=as.double(metastized=="yes"))
df <- df %>% arrange(metastized) # required to align Stan internals and pec 
###############################################################################
sm_exp<- stan_model("~/Documents/Talks/Bayesian_Survival_Models_in_Stan_Meetup_20181105/src/exponential_survival_ppc.stan")
#sm_weibull <- stan_model("~/Documents/Talks/Bayesian_Survival_Models_in_Stan_Meetup_20181105/src/weibull_survival_ppc.stan")
sm_gamma <- stan_model("~/Documents/Talks/Bayesian_Survival_Models_in_Stan_Meetup_20181105/src/gamma_survival_ppc.stan")
#sm_mspline_rw <- stan_model("~/Documents/Talks/Bayesian_Survival_Models_in_Stan_Meetup_20181105/src/survival_parametric_baseline_hazard_rw.stan")
sm_mspline_simplex <- stan_model("~/Documents/Talks/Bayesian_Survival_Models_in_Stan_Meetup_20181105/src/survival_parametric_baseline_hazard_simplex.stan")
###############################################################################
N <- nrow(df)
X <- as.matrix(pull(df, metastized))
is_censored <- pull(df,event)==0
times <- pull(df,time)
time_range <- range(times)
time_min <- time_range[1]
time_max <- time_range[2]
msk_censored <- is_censored == 1
N_censored <- sum(msk_censored)
unique_times <- sort(unique(c(0,times)))
###############################################################################
#ninterior_knots <- 3
#knots <- quantile(times[!msk_censored],head(tail(seq(0,1, length.out = ninterior_knots+2),-1),-1))
#knots <- quantile(times[!msk_censored], probs = c(.1, .5, .9))
knots <- quantile(times[!msk_censored], probs = c(.05, .35, .65, .95))
# We follow Harell's recipe and replace the outer quantiles with the 5 smallest and 5 largest time, respectively
nknots <- length(knots)
knots[1] <- sort(times[!msk_censored])[5]
knots[nknots] <- tail(times[!msk_censored],5)[1]
knots
mspline_degree<- 3
i_spline_basis_evals <- iSpline(times, knots=knots, degree=mspline_degree,
                                intercept=FALSE,Boundary.knots = c(0, max(time_max)))
m_spline_basis_evals <- deriv(i_spline_basis_evals)
i_spline_basis_evals_censored <- i_spline_basis_evals[msk_censored,]
i_spline_basis_evals_uncensored <- i_spline_basis_evals[!msk_censored,]
m_spline_basis_evals_uncensored <- m_spline_basis_evals[!msk_censored,]
nbasis <- dim(i_spline_basis_evals_censored)[2]
###############################################################################
stan_data <- list(N_uncensored=N-N_censored, 
                  N_censored=N_censored, 
                  X_censored=as.matrix(X[msk_censored,]),
                  X_uncensored=as.matrix(X[!msk_censored,]),
                  times_censored=times[msk_censored],
                  times_uncensored = times[!msk_censored],
                  NC=ncol(X),
                  N_times_eval_pec = length(unique_times),
                  times_eval_pec = unique_times,
                  condition=1,
                  m=nbasis,
                  m_spline_basis_evals_uncensored=m_spline_basis_evals_uncensored, 
                  i_spline_basis_evals_uncensored=i_spline_basis_evals_uncensored,
                  i_spline_basis_evals_censored=i_spline_basis_evals_censored
)
###############################################################################
fit_coxph <- coxph(Surv(time, event)~ metastized, data=df, x=TRUE, y=TRUE)
fit_exp <- sampling(sm_exp, data=stan_data, seed=42, chains=4, cores=2, iter=4000)
#fit_weibull <- sampling(sm_weibull, data=stan_data, seed=42, chains=4, cores=2, iter=4000)
fit_gamma <- sampling(sm_gamma, data=stan_data, seed=42, chains=4, cores=2, iter=4000)
#fit_mspline_rw <- sampling(sm_mspline_rw, data=stan_data, seed=42, chains=4, cores=2, iter=4000, control=list(max_treedepth=15, adapt_delta=.99))
fit_mspline_simplex <- sampling(sm_mspline_simplex, data=stan_data, seed=42, chains=4, cores=2, iter=4000)
###############################################################################
post_exp <- as.array(fit_exp)
#post_weibull<- as.array(fit_weibull)
post_gamma <- as.array(fit_gamma)
#post_mspline_rw <- as.array(fit_mspline_rw)
post_mspline_simplex <- as.array(fit_mspline_simplex)
###############################################################################
print(fit_exp, pars=c("betas[1]"))
#print(fit_weibull, pars=c("betas[1]"))
print(fit_gamma, pars=c("betas[1]"))
#print(fit_mspline_rw, pars=c("betas[1]"))
print(fit_mspline_simplex, pars=c("betas[1]"))
###############################################################################
#bayesplot::mcmc_intervals(post_mspline_rw, pars=sprintf("gammas[%d]", 1:nbasis))+xlim(c(0, 0.6))+vline_0()
bayesplot::mcmc_intervals(post_mspline_simplex, pars=sprintf("gammas[%d]", 1:nbasis))+xlim(c(0, 0.6))+vline_0()
###############################################################################
# get confidence intervals for the brier score from posterior samples
brier_scores_exp_df<- get_brier_score_df(post_exp, get_briers_exp_) %>% mutate(model="Exponential")
#brier_scores_weibull_df<- get_brier_score_df(post_weibull, get_briers_weibull_) %>% mutate(model="WeibullBayes")
brier_scores_gamma_df<- get_brier_score_df(post_gamma, get_briers_gamma_) %>% mutate(model="Gamma")

#brier_scores_mspline_rw_df<- get_brier_score_df(post_mspline_rw, get_briers_mspline_rw_) %>% mutate(model="MSpline_RWBayes")
brier_scores_mspline_simplex_df<- get_brier_score_df(post_mspline_simplex, get_briers_mspline_simplex_)%>% mutate(model="MSpline")

brier_scores_coxph_df <- tibble(mean=pec(object=list("model"=fit_coxph),
    formula=Surv(time, event)~1, 
    data=df, exact=TRUE, cens.model="marginal", splitMethod="none",B=300,verbose=FALSE)$AppErr$model,time=unique_times, model="CoxPh",
    low=NA, up=NA
    )

brier_scores_ref_df <- tibble(mean=pec(object=list("model"=fit_coxph),
                                         formula=Surv(time, event)~1, 
                                         data=df, exact=TRUE, cens.model="marginal", splitMethod="none",B=300,verbose=FALSE)$AppErr$Reference,time=unique_times, model="Reference",
                                low=NA, up=NA
)


brier_scores_df <- bind_rows(brier_scores_exp_df, 
                             #brier_scores_weibull_df, 
                             #brier_scores_mspline_rw_df,
                             brier_scores_gamma_df,
                             brier_scores_mspline_simplex_df,
                             brier_scores_coxph_df,
                             brier_scores_ref_df)
ggplot(data=brier_scores_df)+
  geom_hline(yintercept=.25, color='gray', linetype='dashed')+
  geom_stepribbon(mapping=aes(x=time, ymin=low, ymax=up, fill=model),alpha=.1)+
  geom_step(mapping=aes(x=time, y=mean, color=model))+
  ylab("Brier score")+
  xlab("Time")+xlim(c(0, 200))+ylim(c(0, 0.3))+
  scale_fill_manual(values=cbbPalette)+
  scale_colour_manual(values=cbbPalette)
###############################################################################
###############################################################################
# Survival curves
nsamps_per_chain <- dim(post_mspline_simplex)[1]
nchains <- dim(post_mspline_simplex)[2]
times_plot <- 0:time_max
df_surv_mspline_simplex <- purrr::map2_dfr(rep(1:nsamps_per_chain, nchains),rep(1:nchains, each=nsamps_per_chain),
                ~mutate(as.tibble(t(surv_mspline_simplex(times_plot,
                                                  as.matrix(c(0,1)), 
                                                  post_mspline_simplex[.x,.y,sprintf("betas[%d]", 1:ncol(X))],
                                                  post_mspline_simplex[.x,.y,"intercept"],
                                                  post_mspline_simplex[.x,.y,sprintf("gammas[%d]", 1:nbasis)],
                                                  mspline_degree,
                                                  knots,
                                                  attr(i_spline_basis_evals, "Boundary.knots")
                                                  ))
                       ), time=times_plot)
               ) %>%
  group_by(time) %>%
  summarise(mean0=mean(V1) , low0=quantile(V1, probs = c(0.05)), up0=quantile(V1, probs = c(0.975)),
            mean1=mean(V2) , low1=quantile(V2, probs = c(0.05)), up1=quantile(V2, probs = c(0.975))
            ) %>% mutate(model="MSpline")

df_surv_mspline_simplex <- bind_rows(
  mutate(rename(select(df_surv_mspline_simplex, time, mean0, low0, up0, model), mean=mean0, low=low0, up=up0 ), metastized=FALSE),
  mutate(rename(select(df_surv_mspline_simplex, time, mean1, low1, up1, model), mean=mean1, low=low1, up=up1 ), metastized=TRUE)
)
###############################################################################
nsamps_per_chain <- dim(post_exp)[1]
nchains <- dim(post_exp)[2]
times_plot <- 0:time_max
df_surv_exp<- purrr::map2_dfr(rep(1:nsamps_per_chain, nchains),rep(1:nchains, each=nsamps_per_chain),
                                           ~mutate(as.tibble(t(surv_exp(times_plot,
                                                                        as.matrix(c(0,1)), 
                                                                        post_exp[.x,.y,sprintf("betas[%d]", 1:ncol(X))],
                                                                        post_exp[.x,.y,"intercept"]
                                           ))
                                           ), time=times_plot)
) %>%
  group_by(time) %>%
  summarise(mean0=mean(V1) , low0=quantile(V1, probs = c(0.05)), up0=quantile(V1, probs = c(0.975)),
            mean1=mean(V2) , low1=quantile(V2, probs = c(0.05)), up1=quantile(V2, probs = c(0.975))
  ) %>% mutate(model="Exponential")

df_surv_exp <- bind_rows(
  mutate(rename(select(df_surv_exp, time, mean0, low0, up0, model), mean=mean0, low=low0, up=up0 ), metastized=FALSE),
  mutate(rename(select(df_surv_exp, time, mean1, low1, up1, model), mean=mean1, low=low1, up=up1 ), metastized=TRUE)
)

###############################################################################
nsamps_per_chain <- dim(post_gamma)[1]
nchains <- dim(post_gamma)[2]
times_plot <- 0:time_max
df_surv_gamma <- purrr::map2_dfr(rep(1:nsamps_per_chain, nchains),rep(1:nchains, each=nsamps_per_chain),
                              ~mutate(as.tibble(t(surv_gamma(times_plot,
                                                           as.matrix(c(0,1)), 
                                                           post_gamma[.x,.y,sprintf("betas[%d]", 1:ncol(X))],
                                                           post_gamma[.x,.y,"intercept"],
                                                           post_gamma[.x,.y, "alpha"]
                              ))
                              ), time=times_plot)
) %>%
  group_by(time) %>%
  summarise(mean0=mean(V1) , low0=quantile(V1, probs = c(0.05)), up0=quantile(V1, probs = c(0.975)),
            mean1=mean(V2) , low1=quantile(V2, probs = c(0.05)), up1=quantile(V2, probs = c(0.975))
  ) %>% mutate(model="Gamma")

df_surv_gamma <- bind_rows(
  mutate(rename(select(df_surv_gamma, time, mean0, low0, up0, model), mean=mean0, low=low0, up=up0 ), metastized=FALSE),
  mutate(rename(select(df_surv_gamma, time, mean1, low1, up1, model), mean=mean1, low=low1, up=up1 ), metastized=TRUE)
)
###############################################################################
# Frequentist
df_0 <- filter(df, metastized==0)
df_1 <- filter(df_0, event==TRUE)
baseline_hazard_mle <- nrow(df_1)/sum(pull(df_0, "time"))
df_surv_coxph <- tibble(time=times_plot, mean=exp(-times_plot * baseline_hazard_mle), low=NA, up=NA, model="coxph", metastized=FALSE)
df_0 <- filter(df, metastized==1)
df_1 <- filter(df_0, event==TRUE)
baseline_hazard_mle <- nrow(df_1)/sum(pull(df_0, "time"))
df_surv_coxph <- bind_rows(df_surv_coxph, 
                           tibble(time=times_plot, mean=exp(-times_plot * baseline_hazard_mle), low=NA, up=NA, model="coxph", metastized=TRUE)
)

###############################################################################
df_survs <- bind_rows(df_surv_mspline_simplex,df_surv_exp, df_surv_gamma,df_surv_coxph)
ggplot(data=df_survs) +
  geom_ribbon(aes(x=time,ymin = low, ymax = up, fill=model),alpha=.5)+
  geom_line(mapping=aes(x=time, y=mean, color=model, linetype=metastized))+  
  scale_fill_manual(values=cbbPalette)+
  scale_colour_manual(values=cbbPalette)+
  coord_trans(y = "log")
###############################################################################

models <- c("coxph","Exponential","Gamma","MSpline")

get_plot <- function(m) {
  ggplot(data=filter(df_survs, model==m))+
  geom_ribbon(aes(x=time,ymin = low, ymax = up, fill=metastized),alpha=.3)+
  geom_line(mapping=aes(x=time, y=mean,color=metastized, linetype=metastized))+  
  scale_fill_manual(values=cbbPalette)+
  scale_colour_manual(values=cbbPalette)+
  ggtitle(m)+
  coord_trans(y = "log")+
  ylim(c(0.2, 1))
}
plots <- map(models, get_plot)

cowplot::plot_grid(plotlist = plots, align = "h")
###############################################################################
###############################################################################
# PPC's
surv_times_train <- times[!msk_censored]
cens_times_train <- times[msk_censored]
surv_times_rep_exp <- as.matrix(map_dfr(1:dim(post_exp)[2], ~as.tibble(post_exp[,.,sprintf("times_uncensored_sampled[%d]", 1:stan_data$N_uncensored)])))
surv_times_rep_gamma <- as.matrix(map_dfr(1:dim(post_gamma)[2], ~as.tibble(post_gamma[,.,sprintf("times_uncensored_sampled[%d]", 1:stan_data$N_uncensored)])))
# the function picks a random subset of size 500 from the posterior time series (over different chains)
dump_file <- "~/Documents/Talks/Bayesian_Survival_Models_in_Stan_Meetup_20181105/surv_times_rep_mspline_simplex.Rds"
if(!file.exists(dump_file)){
  surv_times_rep_mspline_simplex <- pp_rep_samps_mspline_simplex(1000, as.matrix(X[!msk_censored,]),post_mspline_simplex, mspline_degree, knots, attr(i_spline_basis_evals, "Boundary.knots"), 
                                             nbasis, time_max)
  saveRDS(surv_times_rep_mspline_simplex,dump_file) 
} else {
  surv_times_rep_mspline_simplex <- readRDS(dump_file)
}
################################################################################
p3 <- bayesplot::ppc_dens_overlay(surv_times_train, surv_times_rep_mspline_simplex[sample(1:dim(surv_times_rep_mspline_simplex)[1],50),])+vline_at(knots)+xlim(range(times))+ggtitle("MSpline")
p1 <- bayesplot::ppc_dens_overlay(surv_times_train, surv_times_rep_exp[sample(1:dim(surv_times_rep_exp)[1],50),])+xlim(range(times))+ggtitle("Exponential")
p2 <- bayesplot::ppc_dens_overlay(surv_times_train, surv_times_rep_gamma[sample(1:dim(surv_times_rep_gamma)[1],50),])+xlim(range(times))+ggtitle("Gamma")
p4 <- ggplot(mutate(filter(df,event==TRUE),metastized=as.logical(metastized)), aes(x=time,color=metastized, fill=metastized)) +
      geom_histogram(alpha=0.5, position="identity", binwidth=5)+ggtitle("Histogram")+
      geom_vline(xintercept=mean(pull(filter(df,event==TRUE), "time")),linetype="dotted")
cowplot::plot_grid(p1,p2,p3,p4,align = "h")
################################################################################
color_scheme_set("pink")
bayesplot::bayesplot_grid(
  plots = list(
    ppc_stat(surv_times_train, surv_times_rep_exp, binwidth = 1, stat = "mean"),
    ppc_stat(surv_times_train, surv_times_rep_gamma, binwidth = 1, stat = "mean"),
    ppc_stat(surv_times_train, surv_times_rep_mspline_simplex , binwidth = 1, stat = "mean")
    ),
  titles = c("Exponential","Gamma", "MSpline"),
  legends = FALSE,
  grid_args = list(ncol = 1)
)
################################################################################
color_scheme_set("pink")
bayesplot::bayesplot_grid(
  plots = list(
    ppc_stat(surv_times_train, surv_times_rep_exp, binwidth = 1, stat = "median"),
    ppc_stat(surv_times_train, surv_times_rep_gamma, binwidth = 1, stat = "median"),
    ppc_stat(surv_times_train, surv_times_rep_mspline_simplex , binwidth = 1, stat = "median")
  ),
  titles = c("Exponential","Gamma", "MSpline"),
  legends = FALSE,
  grid_args = list(ncol = 1)
)
################################################################################
color_scheme_set("pink")
bayesplot::bayesplot_grid(
  plots = list(
    ppc_stat(surv_times_train, surv_times_rep_exp, binwidth = 1, stat = "sd"),
    ppc_stat(surv_times_train, surv_times_rep_gamma, binwidth = 1, stat = "sd"),
    ppc_stat(surv_times_train, surv_times_rep_mspline_simplex , binwidth = 1, stat = "sd")
  ),
  titles = c("Exponential","Gamma", "MSpline"),
  legends = FALSE,
  grid_args = list(ncol = 1)
)
################################################################################
color_scheme_set("pink")
bayesplot::bayesplot_grid(
  plots = list(
    ppc_stat(surv_times_train, surv_times_rep_exp, binwidth = 1, stat = "min"),
    ppc_stat(surv_times_train, surv_times_rep_gamma, binwidth = 1, stat = "min"),
    ppc_stat(surv_times_train, surv_times_rep_mspline_simplex , binwidth = 1, stat = "min")
  ),
  titles = c("Exponential","Gamma", "MSpline"),
  legends = FALSE,
  grid_args = list(ncol = 1)
)
################################################################################
color_scheme_set("pink")
bayesplot::bayesplot_grid(
  plots = list(
    ppc_stat(surv_times_train, surv_times_rep_exp, binwidth = 1, stat = "max"),
    ppc_stat(surv_times_train, surv_times_rep_gamma, binwidth = 1, stat = "max"),
    ppc_stat(surv_times_train, surv_times_rep_mspline_simplex , binwidth = 1, stat = "max")
  ),
  titles = c("Exponential","Gamma", "MSpline"),
  legends = FALSE,
  grid_args = list(ncol = 1)
)
################################################################################
km70 <- function(st) {
  tt <- c(st, cens_times_train)
  stat <- c(rep(1, length(st)),rep(0, length(cens_times_train)))
  S <- Surv(tt,stat )
  df_rslt <- tidy(npsurv(S~1))
  filter(df_rslt,time==70)$estimate[1]
}
color_scheme_set("pink")
bayesplot::bayesplot_grid(
  plots = list(
    ppc_stat(surv_times_train, surv_times_rep_exp, stat = "km70"),
    ppc_stat(surv_times_train, surv_times_rep_gamma, stat = "km70"),
    ppc_stat(surv_times_train, surv_times_rep_mspline_simplex , stat = "km70")
  ),
  titles = c("Exponential","Gamma", "MSpline"),
  legends = FALSE,
  grid_args = list(ncol = 1)
)
################################################################################