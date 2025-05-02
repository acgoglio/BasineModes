import numpy as np
import pandas as pd
import xarray as xr
import os

# Load dataset
indir = "/work/cmcc/ag15419/basin_modes/"
ds = xr.open_dataset(os.path.join(indir, "basin_modes_pow_med.nc"))

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

# Rount periods to the second digit
rounded_periods = pd.Series(np.round(periods_in_hours, 2), name="Period")

# Save full list (ordered by frequency)
df_all = rounded_periods.value_counts().reset_index()
df_all.columns = ["Period", "Count"]
df_all["%"] = (df_all["Count"] / df_all["Count"].sum() * 100).round(2)
df_all = df_all.sort_values("Count", ascending=False).reset_index(drop=True)
df_all.to_csv(os.path.join(indir, "periods_all_pow.csv"), index=False)

print("Totale valori prima del filtro:", len(flattened_data))
print("Dopo rimozione NaN:", periods_in_hours.shape[0])
print("Valori min/max:", periods_in_hours.min(), periods_in_hours.max())

# --- Group within ±0.5h ---
sorted_periods = np.sort(periods_in_hours.values)

groups = []
current_group = [sorted_periods[0]]

for p in sorted_periods[1:]:
    if abs(p - current_group[-1]) <= 0.5:
        current_group.append(p)
    else:
        groups.append(current_group)
        current_group = [p]
groups.append(current_group)  # Add the last group

# Summarize groups
group_summary = []
for g in groups:
    center = round(np.mean(g), 2)
    count = len(g)
    group_summary.append((center, count))

df_grouped = pd.DataFrame(group_summary, columns=["Grouped_Period", "Count"])
df_grouped["%"] = (df_grouped["Count"] / df_grouped["Count"].sum() * 100).round(2)
df_grouped = df_grouped.sort_values("Count", ascending=False).reset_index(drop=True)
df_grouped.to_csv(os.path.join(indir, "periods_grouped30min_pow.csv"), index=False)
