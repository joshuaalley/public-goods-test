// Joshua Alley
// Texas A&M University
// Multiple Membership Multilevel Model 


data {
  int<lower = 1> N; // Number of observations
  int<lower = 1> S; // number of states
  int<lower = 1> T; // number of years
  int<lower = 1> A; // number of alliances
  int<lower = 1> M; // number of state-level variables
  int<lower = 1, upper = S> state[N]; // state indicator
  int<lower = 1, upper = T> year[N]; // year indicator
  matrix[N, M] X; // matrix of state-level variables
  matrix[N, A] Z; // matrix of state membership in alliances
  vector[N] y; // outcome 
}

transformed data{
  
  // This section decomposes the sparse matrix Z into a more efficient representation.
  vector[rows(csr_extract_w(Z))] w;
  int v[size(csr_extract_v(Z))]; 
  int u[size(csr_extract_u(Z))]; 
  
  w = csr_extract_w(Z);
  v = csr_extract_v(Z);
  u = csr_extract_u(Z); 
}

parameters {
  real alpha; // overall intercept
  real<lower = 0> sigma; // variance of outcome
  vector[S] alpha_state_std; // better behaved distribution of state intercepts
  vector[T] alpha_year_std; // better behaved distribution of year intercepts
  real<lower = 0> sigma_state; // variance hyperparameter of the state intercepts
  real<lower = 0> sigma_year; // variance hyperparameter of the year intercepts
  vector[M] beta; // vector of state-level coefficients 
  vector[A] lambda; // alliance parameters
  real<lower = 4> nu; // degrees of freedom in t-distribution of outcome

}

transformed parameters {
  vector[S] alpha_state; // state intercepts
  vector[T] alpha_year; // year intercepts
  vector[N] y_hat; // linear prediction of the outcome mean


 alpha_state = 0 + sigma_state * alpha_state_std; // non-centered parameterization, where alpha_state ~ N(0, sigma_state)

alpha_year = 0 + sigma_year * alpha_year_std; // non-centered parameterization, where alpha_state ~ N(0, sigma_state)


// Linear prediction of the state-year spending. csr_matrix_times vector will
// produce a vector as a it multiplies the membership matrix by the vector of alliance characteristics lambda
    y_hat = alpha + alpha_state[state] + alpha_year[year] + csr_matrix_times_vector(N, A, w, v, u, lambda) + X * beta;
    
}

model {
  
  
  alpha ~ normal(0, 3);
  sigma ~ normal(0, 1);
  alpha_year_std ~ normal(0, 1);
  alpha_state_std ~ normal(0, 1); 
  lambda ~ normal(0, 1);
  sigma_state ~ normal(0, 1);
  sigma_year ~ normal(0, 1); 
  beta ~  normal(0, 1);
  nu ~ gamma(2, 0.1); // Prior for degrees of freedom in t-dist
  
  y ~ student_t(nu, y_hat, sigma);
}

generated quantities {
 vector[N] y_pred; //  posterior predictive distribution

 for(i in 1:N)
 y_pred[i] = student_t_rng(nu, y_hat[i], sigma);

}