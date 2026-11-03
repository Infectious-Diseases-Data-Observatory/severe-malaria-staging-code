functions {

  // Function 1: Calculate predicted peripheral stage probabilities from the
  // wrapped Student-t parasite-age distribution under the fixed
  // stage-specific sequestration assumption
  vector wrapped_student_t_probs(real mu_age, real sigma_age, array[] real aa,
                                 int K_stages, vector sequestration_profile) {

    // Store the latent probability assigned to each stage
    vector[K_stages] stage_props;

    for (j in 1:K_stages) {
      stage_props[j] = 0;

      // Sum latent Student-t probability across neighbouring 48-hour cycles
      for (m in -1:1) {
        real upper = (aa[j + 1] + m*48 - mu_age) / sigma_age;
        real lower = (aa[j] + m*48 - mu_age) / sigma_age;

        stage_props[j] += student_t_cdf(upper, 10, 0, 1) -
                          student_t_cdf(lower, 10, 0, 1);
      }

      // Apply the fixed stage-specific sequestration profile
      stage_props[j] *= sequestration_profile[j];
    }

    // Normalise to obtain predicted peripheral stage probabilities
    return stage_props / sum(stage_props);
  }

  // Function 2: Calculate the log-likelihood of the observed peripheral stage
  // counts given the predicted peripheral stage probabilities
  real wrapped_student_t_stage_mix_lpmf(array[] int counts,
                                        real mu1, real sigma1,
                                        real mu2, real sigma2,
                                        real weight1,
                                        array[] real aa,
                                        int K_stages,
                                        vector sequestration_profile) {

    vector[K_stages] p1;
    vector[K_stages] p2;
    vector[K_stages] p;

    // Calculate predicted peripheral stage probabilities for both components
    p1 = wrapped_student_t_probs(
      mu1, sigma1, aa, K_stages, sequestration_profile
    );

    p2 = wrapped_student_t_probs(
      mu2, sigma2, aa, K_stages, sequestration_profile
    );

    // Combine the two predicted peripheral distributions
    p = weight1 * p1 + (1 - weight1) * p2;

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

  // Fixed stage-specific sequestration profile
  vector<lower=0, upper=1>[K_stages] sequestration_profile;

  // Constant describing how "non-centred" the parameterisation of mu_age is
  real<lower=0, upper=1> rho_mu_age;
}

parameters {

  // Mean age is unbounded
  vector[n] mu_age_1;

  // Partially non-centred hierarchy for mean age
  vector[n] z_mu_age;
  real mu_age_mean;
  real<lower=0> mu_age_sd;

  // Patient-specific separation between latent components
  vector<lower=0, upper=24>[n] delta_age;

  // NEW: standard-normal parameters for the non-centred sigma hierarchies
  vector[n] z_log_sigma_centre;
  vector[n] z_log_sigma_difference;

  // Population-level sigma hierarchy
  real log_sigma_centre_mean;
  real<lower=0> log_sigma_centre_sd;

  real log_sigma_difference_mean;
  real<lower=0> log_sigma_difference_sd;

  // Patient-specific mixture weight of mu_age_1
  vector<lower=0, upper=1>[n] weight1;
}

transformed parameters {

  // NEW: reconstruct patient-specific log-sigma parameters
  vector[n] log_sigma_centre;
  vector[n] log_sigma_difference;

  vector[n] log_sigma_age_1;
  vector[n] log_sigma_age_2;

  vector[n] sigma_age_1;
  vector[n] sigma_age_2;

  vector[n] mu_age_2;
  vector[n] mu_age_1_mod48;
  vector[n] mu_age_2_mod48;

  // Non-centred sigma hierarchies
  log_sigma_centre =
    log_sigma_centre_mean +
    log_sigma_centre_sd * z_log_sigma_centre;

  log_sigma_difference =
    log_sigma_difference_mean +
    log_sigma_difference_sd * z_log_sigma_difference;

  // Reconstruct component-specific log sigma
  log_sigma_age_1 =
    log_sigma_centre + 0.5 * log_sigma_difference;

  log_sigma_age_2 =
    log_sigma_centre - 0.5 * log_sigma_difference;

  // Transform back to sigma age
  sigma_age_1 = exp(log_sigma_age_1);
  sigma_age_2 = exp(log_sigma_age_2);

  // Derive second component mean
  mu_age_2 = mu_age_1 + delta_age;

  // Map both means to canonical 0-48 h scale
  for (i in 1:n) {
    mu_age_1_mod48[i] =
      mu_age_1[i] - 48 * floor(mu_age_1[i] / 48);

    mu_age_2_mod48[i] =
      mu_age_2[i] - 48 * floor(mu_age_2[i] / 48);
  }
}

model {

  real combined_sd;

  // Existing partially non-centred mu_age hierarchy
  combined_sd =
    mu_age_sd * sqrt(1 - square(rho_mu_age));

  mu_age_mean ~ normal(12, 24);
  mu_age_sd ~ normal(0, 10);
  z_mu_age ~ normal(0, 1);

  mu_age_1 ~ normal(
    mu_age_mean +
    rho_mu_age * mu_age_sd * z_mu_age,
    combined_sd
  );

  // Patient-specific component separation
  delta_age ~ normal(12, 12);

  // SAME population-level priors as original model
  log_sigma_centre_mean ~ normal(log(4), 1);
  log_sigma_centre_sd ~ normal(0.5, 0.5);

  log_sigma_difference_mean ~ normal(0, 1);
  log_sigma_difference_sd ~ normal(0.5, 0.5);

  // NEW: standard-normal priors generate the same hierarchical
  // distributions through the transformed parameters
  z_log_sigma_centre ~ std_normal();
  z_log_sigma_difference ~ std_normal();

  // Patient-specific mixture weight
  weight1 ~ beta_proportion(0.5, 2);

  for (i in 1:n) {

    target += wrapped_student_t_stage_mix_lpmf(
      stage_obs[i] |
      mu_age_1[i],
      sigma_age_1[i],
      mu_age_2[i],
      sigma_age_2[i],
      weight1[i],
      aa,
      K_stages,
      sequestration_profile
    );
  }
}
