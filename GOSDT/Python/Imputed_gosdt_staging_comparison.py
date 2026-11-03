import pandas as pd

# Paths
WITH_STAGING_FILE = "GOSDT/Python/Imputed_dataset_repeated_validation_by_repeat.csv"
NO_STAGING_FILE = "GOSDT/Python/Imputed_dataset_repeated_validation_by_repeat_no_staging.csv"

# Load results
with_staging = pd.read_csv(WITH_STAGING_FILE)
no_staging = pd.read_csv(NO_STAGING_FILE)

# Check that the same repeats are being compared
print("\nWith-staging rows:", len(with_staging))
print("No-staging rows:", len(no_staging))

print("Same repeats:",
      with_staging["Repeat"].tolist() == no_staging["Repeat"].tolist())
print("Same random states:",
      with_staging["RandomState"].tolist() == no_staging["RandomState"].tolist())

# Merge paired results
comparison = with_staging.merge(no_staging,
                                on = ["Repeat", "RandomState"],
                                suffixes = ("_with_staging", "_no_staging"),
                                validate = "one_to_one")

metrics = ["Accuracy", "BalancedAccuracy", "Sensitivity", "Specificity",
           "Precision", "F1"]

# Calculate paired differences
# Positive differences indicate better performance with staging
# negative differences indicate better performance without staging
for metric in metrics:
    comparison[f"Delta_{metric}"] = (comparison[f"{metric}_with_staging"] -
                                     comparison[f"{metric}_no_staging"])

# Summarise paired differences
summary_rows = []

for metric in metrics:
    delta = comparison[f"Delta_{metric}"]

    summary_rows.append({
        "Metric": metric,
        "Mean_with_staging": comparison[f"{metric}_with_staging"].mean(),
        "Mean_no_staging": comparison[f"{metric}_no_staging"].mean(),
        "Mean_difference": delta.mean(),
        "Median_difference": delta.median(),
        "SD_difference": delta.std(),
        "Minimum_difference": delta.min(),
        "Maximum_difference": delta.max(),
        "Staging_better": int((delta > 0).sum()),
        "Equal": int((delta == 0).sum()),
        "No_staging_better": int((delta < 0).sum())})

summary = pd.DataFrame(summary_rows)

summary["Mean_difference_pp"] = 100 * summary["Mean_difference"]
summary["Median_difference_pp"] = 100 * summary["Median_difference"]
summary["Minimum_difference_pp"] = 100 * summary["Minimum_difference"]
summary["Maximum_difference_pp"] = 100 * summary["Maximum_difference"]

# Save paired results for all 20 repeated train/test partitions
comparison_file = ("GOSDT/Python/Imputed_dataset_staging_comparison_by_repeat.csv")
comparison.to_csv(comparison_file, index = False)

# Save summary of staging versus no-staging performance
summary_file = ("GOSDT/Python/Imputed_dataset_staging_comparison_summary.csv")
summary.to_csv(summary_file, index = False)

print(f"\nSaved '{comparison_file}'")
print(f"Saved '{summary_file}'")

# Display results
print("\nPaired staging versus no-staging comparison\n")
print(summary.to_string(index = False))

print("\nBalanced accuracy differences by repeat "
      "(with staging - without staging)\n")

print(comparison[["Repeat", "RandomState", "BalancedAccuracy_with_staging",
                  "BalancedAccuracy_no_staging",
                  "Delta_BalancedAccuracy"]].to_string(index = False))