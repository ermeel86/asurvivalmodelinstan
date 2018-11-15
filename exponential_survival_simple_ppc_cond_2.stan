functions {
    real exponential_ub_rng(real beta, real ub) {
        real p = exponential_cdf(ub, beta);  // cdf for bounds
        real u = uniform_rng(0, p);
        return (-log1m(u) / beta);  // inverse cdf for value
    }
}
/**************************************************************************************/
data {
    int<lower=0> N_uncensored;                                      
    int<lower=0> N_censored;                                        
    int<lower=1> NC;                                                
    matrix[N_censored,NC] X_censored;                               
    matrix[N_uncensored,NC] X_uncensored;                           
    vector<lower=0>[N_censored] times_censored;                          
    vector<lower=0>[N_uncensored] times_uncensored;                       
}
/**************************************************************************************/
transformed data {
    real max_time;
    real max_time_censored;
    max_time = max(times_uncensored);
    max_time_censored = max(times_censored);
    if(max_time_censored > max_time) max_time = max_time_censored;
}
/**************************************************************************************/
parameters {
    vector[NC] betas;                                     
    real intercept;                                 
}
/**************************************************************************************/
model {
    betas ~ normal(0,2);                                                            
    intercept   ~ normal(-5,2);                                                     
    target += exponential_lpdf(times_uncensored | exp(intercept+X_uncensored*betas)); 
    target += exponential_lccdf(times_censored | exp(intercept+X_censored*betas));  
}
/**************************************************************************************/
generated quantities {
    vector[N_uncensored] times_uncensored_sampled;        
    for(i in 1:N_uncensored) times_uncensored_sampled[i] = exponential_ub_rng(exp(intercept+X_uncensored[i,]*betas), max_time);
}
/**************************************************************************************/
