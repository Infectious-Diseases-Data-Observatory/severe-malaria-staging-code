# HIERARCHICAL MODEL TRACEPLOTS
# COMPLETE CASE

library(rstan)
library(bayesplot)
library(ggplot2)

# ---- Paths ----
PROJECT_DIR <- file.path("Staging/Complete case")
OUTPUT_DIR <- file.path(PROJECT_DIR, "outputs")
RDS_DIR <- file.path(OUTPUT_DIR, "RDS")
TRACEPLOT_DIR <- file.path(OUTPUT_DIR, "traceplots")
dir.create(TRACEPLOT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Input files ----
UNIMODAL_FIT_FILE <- file.path(RDS_DIR, "fit_hierarchical_unimodal.rds")
BIMODAL_FIT_FILE <- file.path(RDS_DIR, "fit_hierarchical_bimodal.rds")

# ---- Load unimodal fit ----
fit_unimodal <- readRDS(UNIMODAL_FIT_FILE)
draws_unimodal <- as.array(fit_unimodal)

# ---- Unimodal population-level diagnostics ----
unimodal_population_pars <- c("mu_age_mean", "mu_age_sd", "log_sigma_age_mean",
                              "log_sigma_age_sd")

unimodal_trace <- bayesplot::mcmc_trace(draws_unimodal, pars = unimodal_population_pars,
                                        facet_args = list(ncol = 1, strip.position = "left"))

unimodal_density <- bayesplot::mcmc_dens_overlay(draws_unimodal,
                                                 pars = unimodal_population_pars)

unimodal_lp_trace <- bayesplot::mcmc_trace(draws_unimodal, pars = "lp__")

ggsave(file.path(TRACEPLOT_DIR, "unimodal_population_traceplot.png"),
       unimodal_trace, width = 8, height = 8, dpi = 300)

ggsave(file.path(TRACEPLOT_DIR, "unimodal_population_density.png"), unimodal_density,
       width = 8, height = 6, dpi = 300)

ggsave(file.path(TRACEPLOT_DIR, "unimodal_log_posterior_traceplot.png"),
       unimodal_lp_trace, width = 8, height = 5, dpi = 300)

rm(fit_unimodal, draws_unimodal)
gc()

# ---- Load bimodal fit ----
fit_bimodal <- readRDS(BIMODAL_FIT_FILE)
draws_bimodal <- as.array(fit_bimodal)

# ---- Bimodal population-level diagnostics ----
bimodal_population_pars_1 <- c("mu_age_mean", "mu_age_sd", "log_sigma_centre_mean")

bimodal_population_pars_2 <- c("log_sigma_centre_sd", "log_sigma_difference_mean",
                               "log_sigma_difference_sd")

bimodal_trace_1 <- bayesplot::mcmc_trace(draws_bimodal, pars = bimodal_population_pars_1,
                                         facet_args = list(ncol = 1, strip.position = "left"))

bimodal_trace_2 <- bayesplot::mcmc_trace(draws_bimodal, pars = bimodal_population_pars_2,
                                         facet_args = list(ncol = 1, strip.position = "left"))

bimodal_density_1 <- bayesplot::mcmc_dens_overlay(draws_bimodal,
                                                  pars = bimodal_population_pars_1)

bimodal_density_2 <- bayesplot::mcmc_dens_overlay(draws_bimodal,
                                                  pars = bimodal_population_pars_2)

bimodal_lp_trace <- bayesplot::mcmc_trace(draws_bimodal, pars = "lp__")

ggsave(file.path(TRACEPLOT_DIR, "bimodal_population_traceplot_1.png"),
       bimodal_trace_1, width = 8, height = 7, dpi = 300)

ggsave(file.path(TRACEPLOT_DIR, "bimodal_population_traceplot_2.png"),
       bimodal_trace_2, width = 8, height = 7, dpi = 300)

ggsave(file.path(TRACEPLOT_DIR, "bimodal_population_density_1.png"),
       bimodal_density_1, width = 8, height = 6, dpi = 300)

ggsave(file.path(TRACEPLOT_DIR, "bimodal_population_density_2.png"),
       bimodal_density_2, width = 8, height = 6, dpi = 300)

ggsave(file.path(TRACEPLOT_DIR, "bimodal_log_posterior_traceplot.png"),
       bimodal_lp_trace, width = 8, height = 5, dpi = 300)

# ---- Bimodal treedepth ----
bimodal_sampler_params <- rstan::get_sampler_params(fit_bimodal, inc_warmup = FALSE)

treedepth_by_chain <- vapply(bimodal_sampler_params,
                             function(x) {sum(x[, "treedepth__"] >= 12)}, numeric(1))

cat("\nTraceplots saved to:\n", TRACEPLOT_DIR, "\n", sep = "")

# ---- Clean up large model objects ----
rm(fit_bimodal, draws_bimodal, bimodal_sampler_params)
gc()