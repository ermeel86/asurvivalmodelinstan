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
    int N = N_uncensored+N_censored;
    max_time = max(times_uncensored);
    max_time_censored = max(times_censored);
    if(max_time_censored > max_time) max_time = max_time_censored;
    
}
/**************************************************************************************/
parameters {
    vector[NC] betas;                                     
    real intercept;
    vector[NC] betas_c;
    real intercept_c;
}
/**************************************************************************************/
model {
    betas ~ normal(0,2);                                                            
    intercept   ~ normal(-5,2); 
    betas_c ~ normal(0,2);
    intercept_c ~ normal(-10,10);

    target += exponential_lpdf(times_uncensored | exp(intercept+X_uncensored*betas)); 
    target += exponential_lccdf(times_censored | exp(intercept+X_censored*betas));  

    target += exponential_lpdf(times_censored | exp(intercept_c+X_censored*betas_c)); 
    target += exponential_lccdf(times_uncensored | exp(intercept_c+X_uncensored*betas_c));  
}
/**************************************************************************************/
generated quantities {
    vector[N] times_samples;
    int<lower=0,upper=1> events[N];
    {
        real samp_e;
        real samp_c;
        for(i in 1:N_uncensored) {
            samp_e = exponential_rng(exp(intercept + X_uncensored[i,]*betas));
            samp_c = exponential_rng(exp(intercept_c + X_uncensored[i,]*betas_c));
            if(samp_e < samp_c) {
                events[i]=1;
                times_samples[i] = samp_e;
            } 
            else {
                events[i] = 0;
                times_samples[i] = samp_c;
            }
        }
        for(i in 1:N_censored) {
            samp_e = exponential_rng(exp(intercept + X_censored[i,]*betas));
            samp_c = exponential_rng(exp(intercept_c + X_censored[i,]*betas_c));
            if(samp_e < samp_c) {
                events[i+N_uncensored]=1;
                times_samples[i+N_uncensored] = samp_e;
            } else {
                events[i+N_uncensored]=0;
                times_samples[i+N_uncensored] = samp_c;
            }
        }
    }
    for(i in 1:N) {
        if(times_samples[i] > max_time) {
            events[i] = 0;
            times_samples[i] = max_time;
        }
    }
    
}
/**************************************************************************************/
