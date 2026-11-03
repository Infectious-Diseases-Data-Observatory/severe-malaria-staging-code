functions {

  // Function 1: Calculate predicted peripheral stage probabilities from the
  // wrapped Student-t parasite-age distribution under the fixed
  // stage-specific sequestration assumption
  vector wrapped_student_t_probs(real mu_age, real sigma_age, array[] real aa, int K_stages,
                                 vector sequestration_profile) {

    // Store the latent probability assigned to each stage
    vector[K_stages] stage_props;

    for (j in 1:K_stages) {
      stage_props[j] = 0;

      // Sum latent Student-t probability across neighbouring 48-hour cycles
      for (m in -1:1) {
        real upper = (aa[j + 1] + m*48 - mu_age) / sigma_age;
        real lower = (aa[j] + m*48 - mu_age) / sigma_age;

        stage_props[j] += student_t_cdf(upper, 10, 0, 1) - student_t_cdf(lower, 10, 0, 1);
      }

      // Apply the fixed stage-specific sequestration profile
      stage_props[j] *= sequestration_profile[j];
    }

    // Normalise to obtain predicted peripheral stage probabilities
    return stage_props / sum(stage_props);
  }

  // Function 2: Calculate the log-likelihood of the observed stage counts given the
  // predicted stage probabilities
  real wrapped_student_t_stage_lpmf(array[] int counts, real mu_age, real sigma_age,
                                    array[] real aa, int K_stages,
                                    vector sequestration_profile) {
    vector[K_stages] p;

    // Calculate predicted peripheral stage probabilities
    p = wrapped_student_t_probs(mu_age, sigma_age, aa, K_stages, sequestration_profile);

    // Likelihood of the observed peripheral stage counts
    return multinomial_lpmf(counts | p);
  }
}

data {
  int<lower=1> K_stages;
  int<lower=1> n;

  // Observed peripheral-blood counts for the seven microscopy stages
  array[n, K_stages] int<lower=0> stage_obs;

  // Parasite-age boundaries defining the seven microscopy stages
  array[K_stages + 1] real aa;

  // Fixed stage-specific sequestration profile.
  // Values represent the relative contribution of each latent parasite stage
  // to the peripheral-blood stage distribution.
  vector<lower=0, upper=1>[K_stages] sequestration_profile;

  // Constant describing how "non-centred" the parameterisation of mu_age is
  // rho_mu_age = 0 is fully centred
  // rho_mu_age = 1 is fully non-centred
  // rho_mu_age = 0.5 is partially non-centred
  real<lower=0, upper=1> rho_mu_age;
}

parameters {

  // Mean age is unbounded
  vector[n] mu_age;

  // Set up hierarchy for mean age
  // z_mu_age is z-score used in the non-centred part of the parameterisation if required
  vector[n] z_mu_age;
  real mu_age_mean;
  real<lower=0> mu_age_sd;

  // Set up hierarchy for sigma age on the log scale
  // Ensures positivity and accounts for right-skewness
  vector[n] log_sigma_age;
  real log_sigma_age_mean;
  real<lower=0> log_sigma_age_sd;
}

transformed parameters {
  vector[n] sigma_age;
  vector[n] mu_age_mod48;

  // Transform log sigma age back to sigma age
  sigma_age = exp(log_sigma_age);

  // Map unbounded mu_age to the canonical 0-48h scale
  for (i in 1:n) {
    mu_age_mod48[i] = mu_age[i] - 48 * floor(mu_age[i] / 48);
  }
}

model {
  real combined_sd;

  // Combine centred and non-centred standard deviations of mu_age
  // Incorporate rho_mu_age
  combined_sd = mu_age_sd * sqrt(1 - square(rho_mu_age));

  // Broad priors for mu_age
  mu_age_mean ~ normal(12, 24);
  mu_age_sd ~ normal(0, 10);
  z_mu_age ~ normal(0, 1);

  // Combine centred and non-centred means of mu_age
  // Incorporate rho_mu_age
  mu_age ~ normal(mu_age_mean + rho_mu_age * mu_age_sd * z_mu_age, combined_sd);

  // Priors centred on typical sigma_age
  log_sigma_age_mean ~ normal(log(4), 1);
  log_sigma_age_sd ~ normal(0.5, 0.5);
  log_sigma_age ~ normal(log_sigma_age_mean, log_sigma_age_sd);

  for (i in 1:n) {

    // Compare the predicted peripheral stage distribution with
    // the observed peripheral stage counts
    target += wrapped_student_t_stage_lpmf(stage_obs[i] | mu_age[i], sigma_age[i], aa,
                                           K_stages, sequestration_profile);
  }
}
