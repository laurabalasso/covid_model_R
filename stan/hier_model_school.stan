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
  
  // apply the above functions and correct the number of positives with exposures 
  matrix corrected_positives(int J, int N, int N_nonzero, int[] nonzero_days, int length_delay, vector delay, matrix conv_gt, matrix r_t, real seed, matrix exposures){
   
    matrix[N_nonzero, J] corrected_positives;
    vector[N] infected ;
    for(j in 1:J){
      infected = get_infected(N, conv_gt, r_t[, j], seed);
      infected = convolve(N, length_delay, infected, delay);
      infected = infected .* exposures[, j];
      corrected_positives[, j] = infected[nonzero_days];
      
    }
  
    return corrected_positives;
  }


}


data {
  int<lower=0> J; // number of groups (regions)
  int<lower=0> N; // total number of observations
  int<lower=0> N_nonzero; // number of observations where total is not 0 
  int<lower=0> nonzero_days[N_nonzero]; // nonzero days index

  matrix[N-1, N] conv_gt; // generation time distribution matrix
  int<lower = 0> length_delay;
  vector[length_delay] p_delay; // delay distribution 

  matrix[N, J] exposures; 
  
  int<lower=0> nonzero_positives[N_nonzero, J]; 
  
  vector<lower=0, upper=1>[N] school ;
  
}


parameters {
  matrix[N, J] log_rt;
  real alpha; // log_rt mean in of groups
  real<lower=0> tau; // between groups variance
  real<lower=0> inv_phi;
  real<lower=0> seed; // initial value for infections vector
  real beta_school;
  
}

transformed parameters{
  real phi = inv(inv_phi);
  matrix<lower = 0>[N, J] r_t ; 
  matrix[N_nonzero, J] eta; // negative binomial location parameters
  

  r_t = exp(log_rt);
  eta = corrected_positives(J,N,N_nonzero, nonzero_days, length_delay, p_delay, conv_gt, r_t, seed, exposures);
}


model {
  
  beta_school ~ normal(1, 10);
  inv_phi ~ normal(0,1);
  seed ~ exponential(1/0.02);
  tau ~ cauchy(0, 2.5);
  
  alpha ~ normal(0, tau);
  
  log_rt[1, ] ~ normal(0, 10);
  for(n in 2:N){
    log_rt[n, ] ~ normal( alpha + log_rt[n-1, ] + beta_school * school[n] , 0.035);
  }
 
  
  for(j in 1:J){
    nonzero_positives[ ,j ] ~ neg_binomial_2(eta[,j], phi);
  }


}


generated quantities {
  
  real log_lik[N_nonzero, J];
  real y_rep[N_nonzero, J];
  
  for(n in 1:N_nonzero){
    for(j in 1:J){
      real eta_nj = eta[n, j];
      y_rep[n, j] = neg_binomial_2_safe_rng(eta_nj, phi);
      log_lik[n, j]  = neg_binomial_2_lpmf(nonzero_positives[n,j]|eta_nj, phi);
      
    }
  }
}


