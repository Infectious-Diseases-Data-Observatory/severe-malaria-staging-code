import pandas as pd

# Read repeated validation results
staging = pd.read_csv("Original_dataset_repeated_validation.csv")
no_staging = pd.read_csv("Original_dataset_repeated_validation_no_staging.csv")

# Keep held-out test results only
staging = staging[staging["Dataset"] == "Test"].copy()
no_staging = no_staging[no_staging["Dataset"] == "Test"].copy()

# Keep the performance measures to compare
metrics = ["Sensitivity", "Specificity", "Accuracy"]

# Merge results from identical train/test splits
comparison = staging[["RandomState"] + metrics].merge(
    no_staging[["RandomState"] + metrics], on = "RandomState",
    suffixes = ("_staging", "_no_staging"), validate = "one_to_one")

# Calculate paired differences:
# positive values mean performance was higher with staging
for metric in metrics:
    comparison[f"{metric}_difference"] = (comparison[f"{metric}_staging"] -
                                          comparison[f"{metric}_no_staging"])

# Save results for every train/test split
comparison.to_csv("Original_dataset_staging_comparison.csv", index = False)

# Summarise paired differences
summary_rows = []

for metric in metrics:
    difference = comparison[f"{metric}_difference"]

    summary_rows.append({"Metric": metric, "MeanDifference": difference.mean(),
                         "SD": difference.std(),
                         "MedianDifference": difference.median(),
                         "MinimumDifference": difference.min(),
                         "MaximumDifference": difference.max(),
                         "StagingBetter": (difference > 0).sum(),
                         "Equal": (difference == 0).sum(),
                         "NoStagingBetter": (difference < 0).sum()})

summary = pd.DataFrame(summary_rows)

# Save summary
summary.to_csv("Original_dataset_staging_comparison_summary.csv", index = False)

print("\nPaired staging versus no-staging comparison:")
print(summary.to_string(index = False))