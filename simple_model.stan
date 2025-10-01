data {
    
  int <lower=0> N; //the number of individuals  
 
  array[N] int <lower=0, upper=1> Y; // Outcome
  real<lower=0> alpha;                  // Beta shape 1
  real<lower=0> beta;                   // Beta shape 2
 
}

parameters {

  real<lower=0, upper=1> annual_foi_raw; // Bernoulli probability

}

   
model {

  annual_foi_raw ~ beta(alpha, beta);

  Y ~ bernoulli(annual_foi_raw);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = bernoulli_lpmf(Y[n] | annual_foi_raw);
  }
}



