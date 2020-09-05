data {
  int<lower=0> TT; // in-sample length
  int<lower=0> oos; // out-of-sample length
  int<lower=1> n; // of space dimensions
  vector[n] y[TT]; // mean corrected returns at time t for n assets
  vector[n] mu; // mean of innovations
  }
parameters {
  cov_matrix[n] H; // covariance of returns
  //matrix[n,n] M; // autoregression for sigma in WAR
  real<lower=0> nu; // degrees of freedom
  //vector[n] mu; // mean of returns in t-distribution
}

model {
  for (t in 1:TT) {
    y[t] ~ multi_student_t(nu, mu, H); 
  }
}
generated quantities {
	vector[n] y_gen[TT];
	matrix[n,n] Sigma[oos];
	 
	for (t in 1:TT) {
	  y_gen[t] = multi_student_t_rng(nu, mu, H); // predict y
	}

	for (t in 1:oos) {
		Sigma[t] = inv_wishart_rng(nu+n+1, H);
	}
}

