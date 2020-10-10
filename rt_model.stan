functions {
  
  // safe Negative Binomial random number generator to avoid overflow problems 
  
  int neg_binomial_2_safe_rng(real eta, real phi) {
    real gamma_rate = gamma_rng(phi, phi / eta);
    if (gamma_rate >= exp(20.79))
      return -9;
      
    return poisson_rng(gamma_rate);
  }
  
  
  // generation of infections vector from r_t and the generation time distribution 
  
  vector get_infected(int N, matrix conv_gt, vector r_t, real seed){
    vector[N] y;
    y[1] = seed;
    for(i in 2:N){
      y[i] = 0;
    }
    
    for(t in 2:N){
      vector[N] gt = to_vector(conv_gt[t-1]);
      y[t] = sum(y .* r_t  .* gt);
    }
    
    return y;
  }
  
  
  vector reverse(int N, vector v){
    vector[N] s ;
    for(i in 1:N){
      s[i] = v[N-i+1];
    }
    return s;
  }
  
  // convolution of infection vector with symptom onset - test delay distribution
  vector convolve(int N, int N_delay, vector infections , vector delay ){
    
    vector[N] out;
    vector[N_delay] rev = reverse(N_delay, delay);
    for(i in 1:N_delay)  {
      out[i] = dot_product(infections[1:i], rev[(N_delay-i+1) : N_delay]);
    }

    for(i in 1:(N-N_delay)){
       out[i+N_delay] = dot_product(infections[i:(i+N_delay-1)], rev);
     }
    
    return out;
    
  }
  
  
  vector corrected_positives(int N, int length_delay, vector delay, matrix conv_gt, vector r_t, real seed, vector exposures, int N_nonzero, int[] nonzero_days){
    vector[N] infected = get_infected(N, conv_gt, r_t, seed);
    vector[N_nonzero] nz_infected;
    infected = convolve(N, length_delay, infected, delay);
    
    nz_infected = (infected .* exposures)[nonzero_days];
    return nz_infected;
  }


}

data {
  int<lower=0> N; // total number of observations
  matrix[N-1, N] conv_gt; // generation time distribution matrix
  int<lower = 0> length_delay;
  vector[length_delay] p_delay; // delay distribution 
  vector[N] exposures; 
  
  int<lower=0> N_nonzero; // number of nonzero total observation
  int<lower=0> nonzero_positives[N_nonzero]; // nonzero vector of positives 
  int<lower=0> nonzero_days[N_nonzero]; // nonzero days index 
  }

parameters {
  vector[N] log_r_t;
  real<lower=0> alpha;
  real<lower=0> seed; // initial value for infections vector
}

transformed parameters{
  vector<lower = 0>[N] r_t = exp(log_r_t);
  vector[N_nonzero] mu = corrected_positives(N, length_delay, p_delay, conv_gt, r_t, seed, exposures, N_nonzero, nonzero_days);

}

model {
  
  alpha ~ gamma(6,1);
  seed ~ exponential(1/0.02);
  log_r_t[1] ~ normal(0, 100) ; 
  
  for(n in 2:N){
    log_r_t[n] ~ normal(log_r_t[n-1], 0.035);
  }
  
  nonzero_positives ~ neg_binomial_2(mu, alpha);
  
}


generated quantities{
 
 
  real y_rep[N_nonzero];
  
  for(n in 1:N_nonzero){
    y_rep[n] = neg_binomial_2_safe_rng(mu[n], alpha);
  }
  
  
}
