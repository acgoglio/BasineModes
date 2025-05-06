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

# --- Greedy grouping algorithm ---
tolerance = 0.4
remaining = rounded_periods.copy()
greedy_groups = []

while not remaining.empty:
    mode = remaining.mode()[0]
    group = remaining[np.abs(remaining - mode) <= tolerance]
    greedy_groups.append((round(group.mean(), 2), len(group)))
    remaining = remaining.drop(group.index)

df_greedy = pd.DataFrame(greedy_groups, columns=["Grouped_Period", "Count"])
df_greedy["%"] = (df_greedy["Count"] / df_greedy["Count"].sum() * 100).round(2)
df_greedy = df_greedy.sort_values("Count", ascending=False).reset_index(drop=True)
df_greedy.to_csv(os.path.join(indir, "periods_grouped_greedy_amp.csv"), index=False)

# --- Histogram of greedy grouped periods ---
plt.figure(figsize=(8, 4))
plt.bar(df_greedy["Grouped_Period"], df_greedy["Count"],
        width=0.4, color="tab:green", edgecolor="black")

plt.xlabel(f"Grouped Period (hours ±{tolerance}h)")
plt.ylabel("Frequency (grid points)")
plt.title("Frequency of Mode Periods in the Mediterranean Sea (Grouped)")
plt.grid(axis="y", linestyle="--", alpha=0.6)
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(os.path.join(indir, "hist_grouped_greedy_amp.png"), dpi=300)
plt.show()
