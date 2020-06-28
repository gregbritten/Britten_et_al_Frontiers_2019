data {
	int<lower=1> N;       // number of observations
	vector[N]    Y;       // response variable
	int<lower=1> K;       // number of population-level effects
	matrix[N, K] X;       // population-level design matrix
	int<lower=1> N_1;     // number of grouping levels
	int<lower=1> M_1;     // number of coefficients per level
	int<lower=1> J_1[N];  // grouping indicator per observation
	vector[N]    Z_1_1;
	vector[N]    Z_1_2;
	int<lower=1> NC_1;    // number of group-level correlations
	int<lower=1> N_2;     // number of grouping levels
	int<lower=1> M_2;     // number of coefficients per level
	int<lower=1> J_2[N];  // grouping indicator per observation
	vector[N]    Z_2_1;
	vector[N]    Z_2_2;
	int<lower=1> NC_2;    // number of group-level correlations
}
transformed data {
	int Kc =      K - 1;
	matrix[N, Kc] Xc;  // centered version of X without an intercept
	vector[Kc]    means_X;  // column means of X before centering
	for (i in 2:K) {
		means_X[i - 1] = mean(X[, i]);
		Xc[, i - 1]    = X[, i] - means_X[i - 1];
	}
}
parameters {
	vector[Kc]                b;         // population-level effects
	real                      Intercept;
	real<lower=0>             sigma;     // residual SD
	vector<lower=0>[M_1]      sd_1;      // group-level standard deviations
	matrix[M_1, N_1]          z_1;       // standardized group-level effects
	cholesky_factor_corr[M_1] L_1;
	vector<lower=0>[M_2]      sd_2;      // group-level standard deviations
	matrix[M_2, N_2]          z_2;       // standardized group-level effects
	cholesky_factor_corr[M_2] L_2;
}
transformed parameters {
	matrix[N_1, M_1] r_1 = (diag_pre_multiply(sd_1, L_1) * z_1)';
	vector[N_1] r_1_1 = r_1[, 1];
	vector[N_1] r_1_2 = r_1[, 2];
	matrix[N_2, M_2] r_2 = (diag_pre_multiply(sd_2, L_2) * z_2)';
	vector[N_2] r_2_1 = r_2[, 1];
	vector[N_2] r_2_2 = r_2[, 2];
}
model {
  vector[N] mu = Intercept + Xc * b;
  for (n in 1:N) {
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n] + r_2_1[J_2[n]] * Z_2_1[n] + r_2_2[J_2[n]] * Z_2_2[n];
  }
  target += uniform_lpdf(b | -3,0);
  target += student_t_lpdf(Intercept | 3, 2, 10);
  target += student_t_lpdf(sigma | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += uniform_lpdf(sd_1[1] | 1E-10,10)
    - 1 * uniform_lccdf(0 | 1E-10,10);
  target += uniform_lpdf(sd_1[2] | 1E-10,3)
    - 1 * uniform_lccdf(0 | 1E-10,3);
  target += normal_lpdf(to_vector(z_1) | 0, 1);
  target += lkj_corr_cholesky_lpdf(L_1 | 1);
  target += uniform_lpdf(sd_2[1] | 1E-10,10)
    - 1 * uniform_lccdf(0 | 1E-10,10);
  target += uniform_lpdf(sd_2[2] | 1E-10,3)
    - 1 * uniform_lccdf(0 | 1E-10,3);
  target += normal_lpdf(to_vector(z_2) | 0, 1);
  target += lkj_corr_cholesky_lpdf(L_2 | 1);

  target += normal_lpdf(Y | mu, sigma);
}
generated quantities {
  real b_Intercept = Intercept - dot_product(means_X, b);
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  corr_matrix[M_2] Cor_2 = multiply_lower_tri_self_transpose(L_2);
  vector<lower=-1,upper=1>[NC_2] cor_2;
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
  for (k in 1:M_2) {
    for (j in 1:(k - 1)) {
      cor_2[choose(k - 1, 2) + j] = Cor_2[j, k];
    }
  }
}
