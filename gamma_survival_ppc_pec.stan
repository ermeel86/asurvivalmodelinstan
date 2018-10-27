/**************************************************************************************/
data {
    int<lower=0> N_uncensored;                                      
    int<lower=0> N_censored;                                        
    int<lower=1> NC;                                                
    int<lower=0> N_times_eval_pec;
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
    real<lower=0> alpha;
}
/**************************************************************************************/
model {
    alpha ~ normal(1, .2);
    betas ~ normal(0,2);                                                            
    intercept   ~ normal(-5,2);                                                     
    target += gamma_lpdf(times_uncensored | alpha, exp(intercept+X_uncensored*betas)); 
    target += gamma_lccdf(times_censored | alpha, exp(intercept+X_censored*betas));  
}
/**************************************************************************************/
generated quantities {
    vector[N_uncensored] times_uncensored_sampled;
    matrix[N,N_times_eval_pec] survs;
    {
        real tmp;
        real max_time;
        real max_time_censored;
        max_time = max(times_uncensored);
        max_time_censored = max(times_censored);
        if(max_time_censored > max_time) max_time = max_time_censored;
        
        for(i in 1:N_uncensored) {
            tmp= max_time + 1; 
            while(tmp > max_time) {
                tmp = gamma_rng(alpha, exp(intercept+X_uncensored[i,]*betas));
            }
            times_uncensored_sampled[i] = tmp;
        }
    }

    for(i in 1:N_uncensored) {
        for(j in 1:N_times_eval_pec) {
            survs[i,j] = 1- gamma_cdf(times_eval_pec[j],alpha, exp(intercept+X_uncensored[i,]*betas));
        }
    }
    for(i in 1:N_censored) {
        for(j in 1:N_times_eval_pec) {
            survs[i+N_uncensored,j] = 1-gamma_cdf(times_eval_pec[j],alpha, exp(intercept+X_censored[i,]*betas));
        }
    }
}
/**************************************************************************************/
