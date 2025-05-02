import numpy as np
import xarray as xr
import pandas as pd
import os

# --- Paths ---
indir = "/work/cmcc/ag15419/basin_modes/"
infile = os.path.join(indir, "basin_modes_amp_med.nc")
csvfile = os.path.join(indir, "periods_grouped30min_amp.csv")
outfile = os.path.join(indir, "mode_groups_amp.nc")

# --- Load NetCDF dataset ---
ds = xr.open_dataset(infile)

# Get list of m?_T variables (mode periods)
modes_vars = [var for var in ds.data_vars if var.startswith("m") and "_T" in var]

# Get lat/lon and shape
nav_lat = ds["nav_lat"].values
nav_lon = ds["nav_lon"].values
field_shape = nav_lat.shape

# --- Load grouped periods from CSV ---
grouped_df = pd.read_csv(csvfile)
group_centers = grouped_df["Grouped_Period"].values

# --- Create output fields ---
fields = {f"mode_{gp:.2f}h": np.zeros(field_shape, dtype=int) for gp in group_centers}

# --- Assign mode numbers ---
for i, var in enumerate(modes_vars):  # i+1 is the mode number
    mode_period_raw = ds[var].values

    # Try to convert mode_period to hours
    if np.issubdtype(mode_period_raw.dtype, np.timedelta64):
        mode_period = mode_period_raw.astype("timedelta64[s]").astype(float) / 3600
    else:
        mode_period = mode_period_raw.astype(float)  # assume already in hours

    # Optional: round to 2 decimals to avoid precision mismatch
    mode_period = np.round(mode_period, 2)

    for gp in group_centers:
        match = np.abs(mode_period - gp) <= 0.5
        if np.any(match):
            print(f"Mode {i+1} matches group {gp:.2f}h in {np.sum(match)} points")
        fields[f"mode_{gp:.2f}h"][match] = i + 1

# --- Save to xarray.Dataset ---
output_ds = xr.Dataset(
    {name: (("y", "x"), data) for name, data in fields.items()},
    coords={"nav_lat": (("y", "x"), nav_lat), "nav_lon": (("y", "x"), nav_lon)}
)

output_ds.to_netcdf(outfile)
print(f"Saved grouped mode maps to: {outfile}")
