/**************************************************************************************/
data {
    int<lower=0> N_uncensored;                                      
    int<lower=0> N_censored;                                        
    int<lower=0> N_times_eval_pec;
    int<lower=1> NC;                                                
    matrix[N_censored,NC] X_censored;                               
    matrix[N_uncensored,NC] X_uncensored;                           
    vector<lower=0>[N_censored] times_censored;                          
    vector<lower=0>[N_uncensored] times_uncensored;
    vector<lower=0>[N_times_eval_pec] times_eval_pec;
}
/**************************************************************************************/
transformed data {
    int<lower=0> N = N_uncensored+N_censored;
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
    target += exponential_lccdf(times_censored  | exp(intercept+X_censored*betas));  
}
/**************************************************************************************/
generated quantities {
    vector[N_uncensored] times_uncensored_sampled;
    matrix[N,N_times_eval_pec] survs;
    for(i in 1:N_uncensored) {
        times_uncensored_sampled[i] = exponential_rng(exp(intercept+X_uncensored[i,]*betas));
    }
    for(i in 1:N_uncensored) {
        for(j in 1:N_times_eval_pec) {
            survs[i,j] = 1- exponential_cdf(times_eval_pec[j],exp(intercept+X_uncensored[i,]*betas));
        }
    }
    for(i in 1:N_censored) {
        for(j in 1:N_times_eval_pec) {
            survs[i+N_uncensored,j] = 1-exponential_cdf(times_eval_pec[j],exp(intercept+X_censored[i,]*betas));
        }
    }
}
/**************************************************************************************/
