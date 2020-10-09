functions {
  
  int neg_binomial_2_safe_rng(real eta, real phi) {
    real gamma_rate = gamma_rng(phi, phi / eta);
    if (gamma_rate >= exp(20.79))
      return -9;
      
    return poisson_rng(gamma_rate);
  }
  
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
  
  vector convolve(int N, int M, vector x, vector y ){
    
    vector[N] out;
    vector[M] rev = reverse(M, y);
    for(i in 1:M)  {
      out[i] = dot_product(x[1:i], rev[(M-i+1) : M]);
    }

    for(i in 1:(N-M)){
       out[i+M] = dot_product(x[i:(i+M-1)], rev);
     }
    
    return out;
    
  }
  
  

}

data {
  int<lower=0> N; // total number of observations
  matrix[N-1, N] conv_gt; // generation time distribution matrix
  int<lower = 0> length_delay;
  vector[length_delay] p_delay; // delay distribution 
  vector[N] exposures; 
  
  int<lower=0> N_nonzero; // number of nonzero observation
  int<lower=0> nonzero_positives[N_nonzero]; // nonzero vector of positives 
  int<lower=0> nonzero_days[N_nonzero]; // nonzero days index 
  }

parameters {
  vector[N] log_r_t;
  real<lower=0> alpha;
  real<lower=0> seed;
}

transformed parameters{
  vector<lower = 0>[N] r_t = exp(log_r_t);
  vector[N] infected = get_infected(N, conv_gt, r_t, seed);
  vector[N] infected_delay = convolve(N, length_delay, infected, p_delay);
  vector[N_nonzero] mu = (infected_delay .* exposures)[nonzero_days];

}

model {
  
  alpha ~ gamma(6,1);
  seed ~ exponential(1/0.02);
  log_r_t[1] ~ normal(0, 0.01) ; 
  
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
