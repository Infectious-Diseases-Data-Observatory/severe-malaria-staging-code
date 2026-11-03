from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# Load imputed data
csv_path = Path("Imputed_dataset_long.csv")
df = pd.read_csv(csv_path, low_memory = False)

# Use one completed dataset for the diagnostic correlation heatmap
df[".imp"] = pd.to_numeric(df[".imp"], errors = "coerce")
df = df[df[".imp"] == 1].copy()

# Variables used in imputation workflow
cols = ["age", "coma", "aconvulsions", "hb", "hct", "aanaemia", "ashockcomp",
        "aacidosis", "temp", "ahypoglyc", "be", "bun", "parasitaemia", "tring",
        "sring", "lring", "etroph", "mtroph", "ltroph", "schizont", "outcome"]

# Trial values to plot separately
trials = ["AQUAMAT", "SEAQUAMAT"]

# First calculate correlations for both trials so they share a scale
corr_matrices = {}
all_offdiag_values = []

for trial in trials:

    # Keep only selected columns and make sure values are numeric where possible
    cols_trial = cols.copy()
    if trial == "SEAQUAMAT":
        cols_trial = [c for c in cols_trial if c != "hrp2_plasma"]

    sub = df.loc[df["trial"] == trial, cols_trial].apply(pd.to_numeric,
                                                         errors = "coerce")

    # Correlation of the actual variables
    corr = sub.corr(method = "pearson")

    # Mask the diagonal so self-correlations do not dominate the figure
    plot_corr = corr.copy()
    np.fill_diagonal(plot_corr.values, np.nan)

    corr_matrices[trial] = plot_corr

    # Use only the visible off-diagonal values to define the scale
    offdiag_values = plot_corr.to_numpy()
    offdiag_values = offdiag_values[~np.isnan(offdiag_values)]
    all_offdiag_values.extend(offdiag_values)

# Common colour scale across both plots, centred on zero
max_abs = np.max(np.abs(all_offdiag_values))
vmin = -max_abs
vmax = max_abs

# Blue -> red colormap
cmap = "RdBu_r"

# Slightly emphasise stronger correlations
norm = TwoSlopeNorm(vmin = vmin, vcenter = 0, vmax = vmax)

# Plot heatmaps
for trial in trials:

    plot_corr = corr_matrices[trial]

    plt.figure(figsize = (10, 8))
    im = plt.imshow(plot_corr, aspect = "auto", cmap = cmap, norm = norm)
    plt.colorbar(im, label = "Correlation")

    plt.xticks(range(len(plot_corr.columns)), plot_corr.columns, rotation = 90)
    plt.yticks(range(len(plot_corr.columns)), plot_corr.columns)
    plt.title(f"Correlation Heatmap - {trial}")
    plt.tight_layout()
    plt.savefig(f"Heatmap_{trial}.png", dpi = 300, bbox_inches = "tight")
    plt.show()