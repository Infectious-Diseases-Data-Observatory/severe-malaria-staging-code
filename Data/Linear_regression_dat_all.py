from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatterMathtext, FuncFormatter, NullFormatter

# Load dataset
csv_path = Path("dat_all.csv")
df = pd.read_csv(csv_path, low_memory = False)

# Remove patients with RBC-total < 93 or > 107
rbc_vars = ["tring", "sring", "lring", "etroph", "mtroph", "ltroph", "schizont"]

df[rbc_vars] = df[rbc_vars].apply(pd.to_numeric, errors = "coerce")
df["rbc_total"] = df[rbc_vars].sum(axis = 1)

df = df[(df["rbc_total"] >= 93) &
        (df["rbc_total"] <= 107)].copy()

df.drop(columns = ["rbc_total"], inplace = True)

# Variables used in plotting
pairs = [("hb", "hct", "Hb vs HCT"),
         ("bun", "be", "BUN vs Base Excess"),
         ("mtroph", "ltroph", "Mature vs Late Trophozoites")]

# Display labels for axes
display_name = {"hb": "Hb [g/dL]", "hct": "HCT [L/L]",
                "bun": "BUN [mmol/L]", "be": "Base Excess [mmol/L]",
                "mtroph": "Mature Trophozoites", "ltroph": "Late Trophozoites"}

trial_colors = {"AQUAMAT": "tab:blue", "SEAQUAMAT": "tab:orange"}

# Label only exact powers of 10 as plain numbers
def plain_decade_label(val, pos = None):
    if val <= 0:
        return ""
    exp = np.log10(val)
    if np.isclose(exp, np.round(exp)):
        return f"{val:g}"
    return ""

# Automatically choose plain decades or powers-of-10 notation
def choose_log_tick_style(x):
    positive = x[x > 0]
    if positive.size == 0:
        return "plain"

    decade_min = np.floor(np.log10(positive.min()))
    decade_max = np.ceil(np.log10(positive.max()))
    span = decade_max - decade_min

    if span <= 4:
        return "plain"
    return "scientific"

# Plot each regression in a separate figure
for xcol, ycol, title in pairs:

    # Keep only selected columns and make sure values are numeric where possible
    sub = df[["Trial", xcol, ycol]].copy()
    sub[[xcol, ycol]] = sub[[xcol, ycol]].apply(pd.to_numeric, errors = "coerce")
    sub = sub.dropna()

    # Plot scatter and regression line
    plt.figure(figsize = (8, 6))

    for trial in ["AQUAMAT", "SEAQUAMAT"]:

        trial_sub = sub[sub["Trial"] == trial]

        # Skip if too few points
        if len(trial_sub) < 2:
            continue

        # Pull out x and y values
        x = trial_sub[xcol].to_numpy()
        y = trial_sub[ycol].to_numpy()

        if xcol == "bun":

            # Keep only positive values for log-scale plots
            mask = x > 0
            x = x[mask]
            y = y[mask]

            if len(x) < 2:
                continue

            # Fit a line in log10 space
            logx = np.log10(x)
            slope, intercept = np.polyfit(logx, y, 1)

            # Line for plotting on log axes
            x_line = np.logspace(np.log10(x.min()), np.log10(x.max()), 200)
            y_line = slope * np.log10(x_line) + intercept

            # Calculate the correlation coefficient and R-squared
            r = np.corrcoef(logx, y)[0, 1]
            r2 = r ** 2

        else:
            # Fit a linear regression line
            slope, intercept = np.polyfit(x, y, 1)
            x_line = np.linspace(x.min(), x.max(), 200)
            y_line = slope * x_line + intercept

            # Calculate the correlation coefficient and R-squared
            r = np.corrcoef(x, y)[0, 1]
            r2 = r ** 2

        # Plot scatter and regression line
        if xcol == "hb" and ycol == "hct":
            plt.scatter(x, y, alpha = 0.12, s = 18, color = trial_colors[trial],
                        label = f"{trial} (n = {len(x)}, R² = {r2:.3f})")
        else:
            plt.scatter(x, y, alpha = 0.35, s = 18, color = trial_colors[trial],
                        label = f"{trial} (n = {len(x)}, R² = {r2:.3f})")

        plt.plot(x_line, y_line, linewidth = 2, color = trial_colors[trial])

    # Add labels and annotation
    plt.xlabel(display_name.get(xcol, xcol))
    plt.ylabel(display_name.get(ycol, ycol))
    plt.title(title)

    ax = plt.gca()

    if xcol == "hb" and ycol == "hct":
        ax.set_xlim(left = 0)
        ax.set_ylim(bottom = 0)
    else:
        ax.margins(x = 0, y = 0)

    if xcol == "bun":

        ax.set_xscale("log")

        # Automatically choose log tick formatting
        tick_style = choose_log_tick_style(sub[xcol].to_numpy())

        ax.xaxis.set_major_locator(LogLocator(base = 10.0))
        ax.xaxis.set_minor_locator(LogLocator(base = 10.0, subs = np.arange(2, 10) * 0.1))
        ax.xaxis.set_minor_formatter(NullFormatter())

        if tick_style == "plain":
            ax.xaxis.set_major_formatter(FuncFormatter(plain_decade_label))
        else:
            ax.xaxis.set_major_formatter(LogFormatterMathtext(base = 10.0))

        # Start the x-axis at 1 mmol/L
        ax.set_xlim(left = 1)

    plt.tight_layout()
    plt.legend()
    plt.savefig(f"LR_{xcol}_vs_{ycol}.png", dpi = 300, bbox_inches = "tight")
    plt.show()