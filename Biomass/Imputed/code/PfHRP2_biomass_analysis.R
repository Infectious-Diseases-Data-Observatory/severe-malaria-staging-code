# PfHRP2 BIOMASS ANALYSIS
# IMPUTED COHORT

library(data.table)
library(ggplot2)

# ---- Shared formatting for supplementary late-stage figures ----
late_stage_theme <- theme_bw(base_size = 12) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11),
        panel.grid.minor = element_blank(), plot.margin = margin(8, 8, 8, 8))

options(scipen = 999)

# ---- Paths ----
BIOMASS_DIR <- file.path("Biomass", "Imputed")
DATA_FILE <- file.path("Data", "dat_all.csv")

STAGE_PROPORTIONS_FILE <- file.path("Staging", "Imputed",
                                    "outputs", "CSV", "stage_proportions.csv")

CSV_DIR <- file.path(BIOMASS_DIR, "CSV")
FIGURE_DIR <- file.path(BIOMASS_DIR, "figures")

dir.create(CSV_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Constants and assumptions ----
# Primary African paediatric PfHRP2 biomass coefficient
TOTAL_BIOMASS_COEFFICIENT <- 7.3e7

# Original Dondorp adult coefficient used for sensitivity analysis
DONDORP_TOTAL_BIOMASS_COEFFICIENT <- 1.2e8

# Blood volume = 80 ml/kg = 80,000 ul/kg
CIRCULATING_BIOMASS_COEFFICIENT <- 80000

OBSERVED_STAGE_NAMES <- c("tring", "sring", "lring", "etroph", "mtroph", "ltroph",
                          "schizont")

OBSERVED_LATE_STAGE_NAMES <- c("etroph", "mtroph", "ltroph", "schizont")

LATENT_STAGE_NAMES <- c("tring_latent", "sring_latent", "lring_latent", "etroph_latent",
                        "mtroph_latent", "ltroph_latent", "schizont_latent")

LATENT_LATE_STAGE_NAMES <- c("etroph_latent", "mtroph_latent", "ltroph_latent",
                             "schizont_latent")

# Assumed stage-specific fractions remaining in peripheral circulation:
# tiny ring, small ring, large ring, early trophozoite, mid trophozoite,
# late trophozoite, schizont
CIRCULATION_PROFILE <- c(1.000, 0.500, 0.250, 0.200, 0.035, 0.027, 0.005)

# ---- Load data ----
dat_all <- fread(DATA_FILE)
stage_proportions <- fread(STAGE_PROPORTIONS_FILE)

# ---- AQUAMAT cohort ----
aquamat <- dat_all[Trial == "AQUAMAT"]

cohort_summary <- aquamat[, .(n = .N, unique_ID = uniqueN(ID),
                              survivors = sum(outcome == 0, na.rm = TRUE),
                              deaths = sum(outcome == 1, na.rm = TRUE))]

cat("\nAQUAMAT cohort\n")
print(cohort_summary)

# ---- Biomass analysis sample ----
biomass <- aquamat[outcome %in% c(0, 1) & is.finite(hrp2_plasma) &
                     is.finite(parasitaemia) & is.finite(hct) & hrp2_plasma > 0 &
                     parasitaemia > 0 & hct > 0 & hct < 100,
                   c("ID", "outcome", "hrp2_plasma", "parasitaemia", "hct",
                     OBSERVED_STAGE_NAMES), with = FALSE]

biomass[, outcome_label := factor(outcome, levels = c(0, 1),
                                  labels = c("Survived", "Died"))]
biomass[, hct_prop := hct / 100]

# ---- Primary paediatric biomass estimates ----
# P_total/kg = 7.3e7 * PfHRP2 * (1 - Hct)
# P_circulating/kg = 80000 * peripheral parasitaemia
# P_sequestered = P_total - P_circulating
biomass[, P_total_per_kg := TOTAL_BIOMASS_COEFFICIENT * hrp2_plasma * (1 - hct_prop)]
biomass[, P_circulating_per_kg := CIRCULATING_BIOMASS_COEFFICIENT * parasitaemia]
biomass[, P_sequestered_per_kg_raw := P_total_per_kg - P_circulating_per_kg]
biomass[, P_sequestered_per_kg := pmax(0, P_sequestered_per_kg_raw)]
biomass[, biomass_ratio := P_total_per_kg / P_circulating_per_kg]

# Retain the raw estimate for analysis; truncate negative estimates to zero
# where a biologically interpretable sequestered fraction is required
biomass[, sequestered_fraction_raw := P_sequestered_per_kg_raw / P_total_per_kg]
biomass[, sequestered_fraction_truncated := pmax(0, sequestered_fraction_raw)]

# ---- Original Dondorp coefficient (sensitivity) ----
biomass[, P_total_Dondorp_per_kg := DONDORP_TOTAL_BIOMASS_COEFFICIENT * hrp2_plasma *
          (1 - hct_prop)]

biomass[, P_sequestered_Dondorp_per_kg := P_total_Dondorp_per_kg - P_circulating_per_kg]
biomass[, biomass_ratio_Dondorp := P_total_Dondorp_per_kg / P_circulating_per_kg]

coefficient_comparison <-
  data.table(formulation = c("African paediatric", "Original Dondorp"),
             coefficient = c(TOTAL_BIOMASS_COEFFICIENT, DONDORP_TOTAL_BIOMASS_COEFFICIENT),
             median_total_per_kg = c(median(biomass$P_total_per_kg),
                                     median(biomass$P_total_Dondorp_per_kg)),
             median_biomass_ratio = c(median(biomass$biomass_ratio),
                                      median(biomass$biomass_ratio_Dondorp)),
             negative_sequestration_n = c(sum(biomass$P_sequestered_per_kg_raw < 0),
                                          sum(biomass$P_sequestered_Dondorp_per_kg < 0)),

  negative_sequestration_percent = c(100 * mean(biomass$P_sequestered_per_kg_raw < 0),
                                     100 * mean(biomass$P_sequestered_Dondorp_per_kg < 0)),

             population_circulating_fraction = c(sum(biomass$P_circulating_per_kg) /
                                                   sum(biomass$P_total_per_kg),
                                                 sum(biomass$P_circulating_per_kg) /
                                                   sum(biomass$P_total_Dondorp_per_kg)))

cat("\nBiomass coefficient sensitivity\n")
print(coefficient_comparison)

# ---- Primary biomass summary ----
biomass_summary <- biomass[, .(n = .N,
                               median_total_per_kg = median(P_total_per_kg),
                               median_circulating_per_kg = median(P_circulating_per_kg),
                               median_sequestered_per_kg = median(P_sequestered_per_kg_raw),
                               median_biomass_ratio = median(biomass_ratio),
                               negative_sequestration_n = sum(P_sequestered_per_kg_raw < 0),
                  negative_sequestration_percent = 100 * mean(P_sequestered_per_kg_raw < 0))]

cat("\nPrimary paediatric biomass summary:\n")
print(biomass_summary)

# ---- Figure 3 adaptation (total biomass distribution) ----
p_total_distribution <- ggplot(biomass, aes(x = P_total_per_kg)) +
  geom_histogram(bins = 40, colour = "black", fill = "grey80") +
  scale_x_log10(labels = scales::label_scientific()) +
  labs(title = "Distribution of PfHRP2-derived total parasite biomass",
       subtitle = paste0("AQUAMAT; paediatric formulation; n = ", nrow(biomass)),
       x = "Total parasite biomass (parasites/kg; log10 scale)",
       y = "Patients") +
  theme_bw()

print(p_total_distribution)

ggsave(file.path(FIGURE_DIR, "Figure3_total_biomass_distribution.png"),
       p_total_distribution, width = 8, height = 6, dpi = 300)

# ---- Figure 4 adaptation (biomass by survival outcome) ----
outcome_counts <- biomass[, .N, by = outcome_label]

outcome_biomass_summary <- biomass[, .(n = .N,
                                       median_total = median(P_total_per_kg),
                                       median_circulating = median(P_circulating_per_kg),
                                  median_sequestered = median(P_sequestered_per_kg_raw)),
                                  by = outcome_label]

total_outcome_test <- wilcox.test(P_total_per_kg ~ outcome_label, data = biomass,
                                  exact = FALSE)

circulating_outcome_test <- wilcox.test(P_circulating_per_kg ~ outcome_label,
                                        data = biomass, exact = FALSE)

sequestered_outcome_test <- wilcox.test(P_sequestered_per_kg_raw ~ outcome_label,
                                        data = biomass, exact = FALSE)

outcome_test_summary <- data.table(biomass = c("Total", "Circulating", "Sequestered"),
                                   median_survived = c(biomass[outcome == 0,
                                                               median(P_total_per_kg)],
                                                       biomass[outcome == 0,
                                                               median(P_circulating_per_kg)],
                                                       biomass[outcome == 0,
                                                       median(P_sequestered_per_kg_raw)]),

                                   median_died = c(biomass[outcome == 1,
                                                           median(P_total_per_kg)],
                                                   biomass[outcome == 1,
                                                           median(P_circulating_per_kg)],
                                                   biomass[outcome == 1,
                                                   median(P_sequestered_per_kg_raw)]),

                                   wilcoxon_p = c(total_outcome_test$p.value,
                                                  circulating_outcome_test$p.value,
                                                  sequestered_outcome_test$p.value))

outcome_test_summary[, died_to_survived_ratio := median_died / median_survived]

cat("\nBiomass by outcome:\n")
print(outcome_counts)
print(outcome_biomass_summary)
print(outcome_test_summary)

biomass_outcome_long <- rbindlist(list(biomass[, .(ID, outcome_label,
                                                   biomass_type = "Total",
                                                   biomass_per_kg = P_total_per_kg)],
                                       biomass[, .(ID, outcome_label,
                                                   biomass_type = "Circulating",
                                                   biomass_per_kg = P_circulating_per_kg)]))

biomass_outcome_long[, biomass_type := factor(biomass_type, levels = c("Circulating",
                                                                       "Total"))]

# ---- Plot estimated biomass by survival outcome ----
p_outcome <- ggplot(biomass_outcome_long,
                    aes(x = outcome_label, y = biomass_per_kg, fill = biomass_type)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = c("Circulating" = "#F8766D", "Total" = "#00BFC4"),
                    name = NULL) +
  scale_y_log10(labels = scales::label_scientific()) +
  labs(x = NULL, y = "Parasite biomass (parasites/kg)") + theme_bw(base_size = 12) +
  theme(axis.title.y = element_text(size = 12), axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11), legend.text = element_text(size = 11),
        legend.position = "top", panel.grid.minor = element_blank(),
        plot.margin = margin(8, 8, 8, 8))

print(p_outcome)

ggsave(file.path(FIGURE_DIR, "Figure4_biomass_by_outcome.png"), plot = p_outcome,
       width = 7, height = 5, units = "in", dpi = 300)

# ---- Observed late-stage proportion ----
biomass[, observed_stage_total := rowSums(.SD), .SDcols = OBSERVED_STAGE_NAMES]

# Late stages = early trophozoites + mid trophozoites + late trophozoites + schizonts
biomass[, observed_late_stage_percent := fifelse(observed_stage_total > 0,
                                                 100 * rowSums(.SD) / observed_stage_total,
                                                 NA_real_),
        .SDcols = OBSERVED_LATE_STAGE_NAMES]

# ---- Figure S1 adaptation (full observed-staging sample) ----
observed_stage_data <- biomass[is.finite(observed_late_stage_percent) &
                                 is.finite(sequestered_fraction_raw)]

observed_cor_raw <- cor.test(observed_stage_data$observed_late_stage_percent,
                             observed_stage_data$sequestered_fraction_raw,
                             method = "spearman", exact = FALSE)

observed_cor_truncated <- cor.test(observed_stage_data$observed_late_stage_percent,
                                   observed_stage_data$sequestered_fraction_truncated,
                                   method = "spearman", exact = FALSE)

observed_positive <- observed_stage_data[biomass_ratio >= 1]

observed_cor_positive <- cor.test(observed_positive$observed_late_stage_percent,
                                  observed_positive$sequestered_fraction_raw,
                                  method = "spearman", exact = FALSE)

observed_sensitivity <- data.table(sample = "Full observed-staging sample",
                                   stage_measure = "Observed peripheral late stages",
                                   sequestration_treatment = c("Raw",
                                                 "Negative estimates set to zero",
                                                 "P_total >= P_circulating only"),

                                   n = c(nrow(observed_stage_data),
                                         nrow(observed_stage_data), nrow(observed_positive)),

                                   spearman_rho = c(unname(observed_cor_raw$estimate),
                                                    unname(observed_cor_truncated$estimate),
                                                    unname(observed_cor_positive$estimate)),

                                   p_value = c(observed_cor_raw$p.value,
                                               observed_cor_truncated$p.value,
                                               observed_cor_positive$p.value))

cat("\nFigure S1 adaptation:\n")
print(observed_sensitivity)

# ---- Figure S3(a): Observed peripheral late stages ----
p_observed_late <- ggplot(observed_stage_data,
                          aes(x = observed_late_stage_percent,
                              y = 100 * sequestered_fraction_truncated)) +
  geom_point(alpha = 0.30, size = 0.8) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE) +
  scale_x_continuous(breaks = seq(0, 100, 25), limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(x = "Observed peripheral late stages (%)",
       y = "PfHRP2-derived sequestered biomass (% of total)") +
  late_stage_theme

print(p_observed_late)

ggsave(file.path(FIGURE_DIR, "FigureS1_observed_late_stages.png"),
       plot = p_observed_late, width = 7, height = 5, units = "in", dpi = 300)

# ---- Match biomass and hierarchical staging ----
matched <- merge(biomass, stage_proportions[trial == "AQUAMAT",
                                            c("ID", "model", "mu_age", "sigma_age",
                                              LATENT_STAGE_NAMES), with = FALSE],
                 by = "ID", all = FALSE)

latent_matrix <- as.matrix(matched[, ..LATENT_STAGE_NAMES]) / 100

matched[, latent_late_stage_percent := rowSums(.SD), .SDcols = LATENT_LATE_STAGE_NAMES]
matched[, circulating_fraction_PfHRP2 := P_circulating_per_kg / P_total_per_kg]

# Predicted circulating fraction = sum of each latent stage proportion
# multiplied by its assumed fraction remaining in peripheral circulation
matched[, circulating_fraction_hierarchical := as.numeric(latent_matrix %*%
                                                            CIRCULATION_PROFILE)]

# ---- Matched observed vs latent comparison ----
# Restrict both analyses to the same patients so that differences between
# observed and latent correlations are not caused by different sample sizes
matched_stage_data <- matched[is.finite(observed_late_stage_percent) &
                                is.finite(latent_late_stage_percent) &
                                is.finite(sequestered_fraction_raw)]

stage_sensitivity <- function(data, stage_variable, stage_label) {
  raw_test <- cor.test(data[[stage_variable]], data$sequestered_fraction_raw,
                       method = "spearman", exact = FALSE)
  truncated_test <- cor.test(data[[stage_variable]], data$sequestered_fraction_truncated,
                             method = "spearman", exact = FALSE)

  positive_data <- data[biomass_ratio >= 1]
  positive_test <- cor.test(positive_data[[stage_variable]],
                            positive_data$sequestered_fraction_raw,
                            method = "spearman", exact = FALSE)

  data.table(sample = "Matched observed/latent sample",
             stage_measure = stage_label,
             sequestration_treatment = c("Raw", "Negative estimates set to zero",
                                         "P_total >= P_circulating only"),
             n = c(nrow(data),  nrow(data), nrow(positive_data)),
             spearman_rho = c(unname(raw_test$estimate), unname(truncated_test$estimate),
                              unname(positive_test$estimate)),
             p_value = c(raw_test$p.value, truncated_test$p.value, positive_test$p.value))
}

matched_observed_sensitivity <- stage_sensitivity(matched_stage_data,
                                                  "observed_late_stage_percent",
                                                  "Observed peripheral late stages")

matched_latent_sensitivity <- stage_sensitivity(matched_stage_data,
                                                "latent_late_stage_percent",
                                                "Hierarchical latent late stages")

stage_comparison <- rbindlist(list(observed_sensitivity, matched_observed_sensitivity,
                                   matched_latent_sensitivity), use.names = TRUE)

cat("\nObserved vs latent stage comparison:\n")
print(stage_comparison)

# ---- Hierarchical extension figure ----
matched_latent_cor_truncated <- cor.test(matched_stage_data$latent_late_stage_percent,
                                         matched_stage_data$sequestered_fraction_truncated,
                                         method = "spearman", exact = FALSE)

# ---- Figure S3(b): Inferred latent late stages ----
p_latent_late <- ggplot(matched_stage_data,
                        aes(x = latent_late_stage_percent,
                            y = 100 * sequestered_fraction_truncated)) +
  geom_point(alpha = 0.30, size = 0.8) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE) +
  scale_x_continuous(breaks = seq(0, 100, 25), limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(x = "Inferred latent late stages (%)",
       y = "PfHRP2-derived sequestered biomass (% of total)") +
  late_stage_theme

print(p_latent_late)

ggsave(file.path(FIGURE_DIR, "Extension_latent_late_stages.png"),
       plot = p_latent_late, width = 7, height = 5, units = "in", dpi = 300)

# ---- PfHRP2 vs hierarchical circulating fraction ----
# Compare the independently PfHRP2-derived circulating fraction with the
# circulating fraction implied by the hierarchical staging framework
population_PfHRP2 <- sum(matched$P_circulating_per_kg) / sum(matched$P_total_per_kg)

population_hierarchical <- weighted.mean(matched$circulating_fraction_hierarchical,
                                         w = matched$P_total_per_kg)

population_comparison <- data.table(method = c("PfHRP2-derived", "Hierarchical staging"),
                                    circulating_fraction = c(population_PfHRP2,
                                                             population_hierarchical))

population_comparison[, circulating_percent := 100 * circulating_fraction]

pearson_circulation <- cor.test(matched$circulating_fraction_PfHRP2,
                                matched$circulating_fraction_hierarchical,
                                method = "pearson")

spearman_circulation <- cor.test(matched$circulating_fraction_PfHRP2,
                                 matched$circulating_fraction_hierarchical,
                                 method = "spearman", exact = FALSE)

circulation_agreement <-
  data.table(n = nrow(matched),
             median_PfHRP2 = median(matched$circulating_fraction_PfHRP2),
             median_hierarchical = median(matched$circulating_fraction_hierarchical),
             population_PfHRP2 = population_PfHRP2,
             population_hierarchical = population_hierarchical,
             hierarchical_to_PfHRP2_ratio = population_hierarchical / population_PfHRP2,
             pearson_r = unname(pearson_circulation$estimate),
             pearson_p = pearson_circulation$p.value,
             spearman_rho = unname(spearman_circulation$estimate),
             spearman_p = spearman_circulation$p.value)

circulation_by_model <- matched[, .(n = .N,
                        median_PfHRP2 = median(circulating_fraction_PfHRP2),
                        median_hierarchical = median(circulating_fraction_hierarchical)),
                        by = model]

cat("\nPfHRP2 vs hierarchical circulating fraction:\n")
print(population_comparison)
print(circulation_agreement)
print(circulation_by_model)

# ---- Matched-sample coefficient sensitivity ----
matched[, circulating_fraction_PfHRP2_Dondorp := P_circulating_per_kg /
          P_total_Dondorp_per_kg]

matched_coefficient_sensitivity <-
  data.table(formulation = c("African paediatric", "Original Dondorp"),

             population_circulating_fraction = c(sum(matched$P_circulating_per_kg) /
                                                   sum(matched$P_total_per_kg),
                                                 sum(matched$P_circulating_per_kg) /
                                                   sum(matched$P_total_Dondorp_per_kg)))

matched_coefficient_sensitivity[,
                    hierarchical_population_fraction := population_hierarchical]

matched_coefficient_sensitivity[,
                    hierarchical_to_PfHRP2_ratio := hierarchical_population_fraction /
                      population_circulating_fraction]

cat("\nMatched-sample coefficient sensitivity:\n")
print(matched_coefficient_sensitivity)

# ---- Export results ----
fwrite(biomass, file.path(CSV_DIR, "PfHRP2_biomass.csv"))
fwrite(biomass_summary, file.path(CSV_DIR, "Biomass_summary.csv"))
fwrite(outcome_test_summary, file.path(CSV_DIR, "Outcome_comparison.csv"))
fwrite(stage_comparison, file.path(CSV_DIR, "Stage_comparison.csv"))

fwrite(population_comparison, file.path(CSV_DIR,
                                        "Population_circulation_comparison.csv"))

fwrite(circulation_agreement, file.path(CSV_DIR, "Circulation_agreement.csv"))
fwrite(circulation_by_model, file.path(CSV_DIR, "Circulation_by_model.csv"))
fwrite(coefficient_comparison, file.path(CSV_DIR, "Coefficient_comparison.csv"))

fwrite(matched_coefficient_sensitivity,
       file.path(CSV_DIR, "Matched_coefficient_sensitivity.csv"))

# ---- Complete ----
cat("\nBiomass analysis complete\n")

cat("\nPrimary formulation: African paediatric coefficient = 7.3\n",
    "Sensitivity formulation: original Dondorp coefficient = 12\n", sep = "")

cat("\nSaved CSV outputs to:\n", CSV_DIR, "\n", sep = "")
cat("\nSaved figures to:\n", FIGURE_DIR, "\n", sep = "")