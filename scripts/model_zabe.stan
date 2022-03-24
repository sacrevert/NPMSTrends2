data{
  int<lower=1> N1;
  int<lower=1> visits1[N1];
  int<lower=1> pos1[N1];
  int<lower=1> n_obs1;
  int<lower=0, upper=1> x[n_obs1];
  
  int<lower=1> N2;
  int<lower=1> visits2[N2];
  int<lower=1> pos2[N2];
  int<lower=1> n_obs2;
  vector[n_obs2] logit_u;
  int<lower=1> N1_ind[N2];
}

parameters{
  real<lower=0, upper=1> psi;
  real<lower=0, upper=1> mu;
  real<lower=0> phi;
  real gamma0;
  real<lower=0> gamma1;
  real<lower=0> sigma;
  real<lower=0, upper=1> y[N1];
}

transformed parameters{
  real alpha;
  real beta;
  real logit_p[N1];
  real phi_t;
  
  alpha = phi*mu;
  beta = phi*(1-mu);
  for(i in 1:N1){
    logit_p[i] = gamma0 + gamma1*y[i];
  }
  phi_t = pow(phi, -0.5);
}

model{
  real log_psi;
  real log1m_psi;
  real logit_y[N1];
  
  log_psi = log(psi);
  log1m_psi = log1m(psi);
  
  y ~ beta(alpha, beta);
  
  for(i in 1:N1){
    int x_temp[visits1[i]];
    real logit_p_temp;
    
    logit_y[i] = logit(y[i]);
    x_temp = segment(x, pos1[i], visits1[i]);
    logit_p_temp = logit_p[i];
    
    if(sum(x_temp)>0)
      target += log_psi +
        bernoulli_logit_lpmf(x_temp | logit_p_temp);
    else
      target += log_sum_exp(log_psi + bernoulli_logit_lpmf(x_temp | logit_p_temp),
      log1m_psi);
  }
  
  // for(i in 1:N2){
  //   vector[visits2[i]] logit_u_temp;
  //   logit_u_temp = segment(logit_u, pos2[i], visits2[i]);
  //   //logit_u_temp ~ normal(logit_y[N1_ind[i]], sigma);
  //   target += normal_lpdf(logit_u_temp | logit_y[N1_ind[i]], sigma);
  // }
  
  for(i in 1:N2){
    for(j in 1:visits2[i]){
      logit_u[(pos2[i] + j - 1)] ~ normal(logit_y[N1_ind[i]], sigma);
    }
  }
  
  psi ~ uniform(0, 1);
  mu ~ uniform(0, 1);
  phi ~ cauchy(0, 10);
  gamma0 ~ cauchy(0, 10);
  gamma1 ~ cauchy(0, 5);
  sigma ~ cauchy(0, 5);
}
