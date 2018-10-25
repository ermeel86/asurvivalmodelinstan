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
parameters {
    vector[NC] betas;                                     
    real intercept;
    real<lower=0> alpha;
}
/**************************************************************************************/
model {
    alpha ~ normal(1, .2);
    betas ~ normal(0,2);                                                            
    intercept   ~ normal(-5,2);                                                     
    target += weibull_lpdf(times_uncensored | alpha, exp(intercept+X_uncensored*betas)); 
    target += weibull_lccdf(times_censored | alpha, exp(intercept+X_censored*betas));  
}
/**************************************************************************************/
generated quantities {
    vector[N_uncensored] times_uncensored_sampled;
    for(i in 1:N_uncensored) {
        times_uncensored_sampled[i] = weibull_rng(alpha, exp(intercept+X_uncensored[i,]*betas));
    }
}
/**************************************************************************************/
