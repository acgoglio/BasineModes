import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
mpl.use("Agg")  # For non-interactive backend

# Load dataset
indir = "/work/cmcc/ag15419/basin_modes/"
ds = xr.open_dataset(os.path.join(indir, "basin_modes_amp_med.nc"))

# Extract m?_T variables
modes_vars = [var for var in ds.data_vars if var.startswith("m") and "_T" in var]
mode_data = np.stack([ds[var].values for var in modes_vars], axis=-1)
flattened_data = mode_data.flatten()

# Convert to hours
periods_in_hours = pd.to_timedelta(flattened_data).total_seconds() / 3600

# Drop NaN and invalid values
periods_in_hours = pd.to_numeric(periods_in_hours, errors="coerce")
periods_in_hours = periods_in_hours.dropna()
periods_in_hours = periods_in_hours[periods_in_hours > 0]
periods_in_hours = periods_in_hours[periods_in_hours <= 40]

if periods_in_hours.empty:
    print("Nessun valore valido di periodo tra 0 e 40 ore trovato.")
    exit()

# Round periods to second digit
rounded_periods = pd.Series(np.round(periods_in_hours, 2), name="Period")

# Save full list (ordered by frequency)
df_all = rounded_periods.value_counts().reset_index()
df_all.columns = ["Period", "Count"]
df_all["%"] = (df_all["Count"] / df_all["Count"].sum() * 100).round(2)
df_all = df_all.sort_values("Count", ascending=False).reset_index(drop=True)
df_all.to_csv(os.path.join(indir, "periods_all_amp.csv"), index=False)

print("Totale valori prima del filtro:", len(flattened_data))
print("Dopo rimozione NaN:", periods_in_hours.shape[0])
print("Valori min/max:", periods_in_hours.min(), periods_in_hours.max())

# --- Histogram of all periods ---
plt.figure(figsize=(8, 4))
bars = plt.bar(df_all["Period"], df_all["Count"],
               width=0.06, color="tab:orange", edgecolor="black")

plt.xlabel("Period (hours)")
plt.ylabel("Frequency (grid points)")
plt.title("Frequency of Mode Periods in the Mediterranean Sea")
plt.grid(axis="y", linestyle="--", alpha=0.6)
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(os.path.join(indir, "hist_all_amp.png"), dpi=300)
plt.show()

# --- Group into ±0.5h bins centered on integer hours ---
bins = np.arange(0.5, 40.5, 1.0)  # Bin edges: 0.5, 1.5, ..., 39.5
labels = np.arange(1, 40)         # Bin centers: 1, 2, ..., 39

binned = pd.cut(periods_in_hours, bins=bins, labels=labels)
df_grouped = binned.value_counts().sort_index().reset_index()
df_grouped.columns = ["Grouped_Period", "Count"]
df_grouped["Grouped_Period"] = df_grouped["Grouped_Period"].astype(float)
df_grouped["%"] = (df_grouped["Count"] / df_grouped["Count"].sum() * 100).round(2)
df_grouped.to_csv(os.path.join(indir, "periods_grouped30min_amp.csv"), index=False)

# --- Histogram of grouped periods ---
plt.figure(figsize=(8, 4))
plt.bar(df_grouped["Grouped_Period"], df_grouped["Count"],
        width=0.9, color="tab:blue", edgecolor="black")

plt.xlabel("Grouped Period (hours ±0.5h)")
plt.ylabel("Frequency (grid points)")
plt.title("Frequency of Mode Periods in the Mediterranean Sea")
plt.grid(axis="y", linestyle="--", alpha=0.6)
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(os.path.join(indir, "hist_grouped30min_amp.png"), dpi=300)
plt.show()

# --- Group into ±1h bins centered on integer hours ---
bins = np.arange(0, 41, 2.0)  # → bin edges: 0, 2, 4, ..., 40
labels = np.arange(1, 40, 2)  # → bin centers: 1, 3, 5, ..., 39

binned = pd.cut(periods_in_hours, bins=bins, labels=labels, right=False)
df_grouped = binned.value_counts().sort_index().reset_index()
df_grouped.columns = ["Grouped_Period", "Count"]
df_grouped["Grouped_Period"] = df_grouped["Grouped_Period"].astype(float)
df_grouped["%"] = (df_grouped["Count"] / df_grouped["Count"].sum() * 100).round(2)
df_grouped.to_csv(os.path.join(indir, "periods_grouped1h_amp.csv"), index=False)

# --- Histogram of grouped periods ---
plt.figure(figsize=(8, 4))
plt.bar(df_grouped["Grouped_Period"], df_grouped["Count"],
        width=0.9, color="tab:blue", edgecolor="black")

plt.xlabel("Grouped Period (hours ±1h)")
plt.ylabel("Frequency (grid points)")
plt.title("Frequency of Mode Periods in the Mediterranean Sea")
plt.grid(axis="y", linestyle="--", alpha=0.6)
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(os.path.join(indir, "hist_grouped1h_amp.png"), dpi=300)
plt.show()
