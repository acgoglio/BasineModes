import numpy as np
import xarray as xr
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# --- Paths ---
indir = "/work/cmcc/ag15419/basin_modes/"
infile = os.path.join(indir, "basin_modes_amp_med.nc")
csvfile = os.path.join(indir, "periods_grouped30min_amp.csv")
outfile = os.path.join(indir, "mode_groups_amp.nc")
mesh_mask_file = "/work/cmcc/ag15419/VAA_paper/DATA0/mesh_mask.nc"
output_plot_dir = os.path.join(indir, "mode_plots")  # Directory for saving plots

# Create output plot directory if it doesn't exist
os.makedirs(output_plot_dir, exist_ok=True)

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

# --- Load land-sea mask (tmask) ---
tmask_ds = xr.open_dataset(mesh_mask_file)
tmask = tmask_ds["tmask"].values  # assuming tmask is in the variable "tmask"

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
            print(f"Mode {i} matches group {gp:.2f}h in {np.sum(match)} points")
        fields[f"mode_{gp:.2f}h"][match] = i  # Mode numbers from 0 to 8

# --- Save to xarray.Dataset ---
output_ds = xr.Dataset(
    {name: (("y", "x"), data) for name, data in fields.items()},
    coords={"nav_lat": (("y", "x"), nav_lat), "nav_lon": (("y", "x"), nav_lon)}
)

output_ds.to_netcdf(outfile)
print(f"Saved grouped mode maps to: {outfile}")

# --- Plotting: Mode Group Maps with Land-Sea Mask Contour ---
for gp in group_centers:
    mode_field = fields[f"mode_{gp:.2f}h"]

    plt.figure(figsize=(10, 8))
    plt.contourf(nav_lon, nav_lat, mode_field, levels=np.arange(0, 9), cmap="Reds_r", extend="both")
    
    # Overlay land-sea mask as black contours
    contour = plt.contour(nav_lon, nav_lat, tmask[0, :, :], levels=[0.5], colors="black", linewidths=1)

    plt.title(f"Mode Groups for Period: {gp:.2f}h")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.colorbar(label="Mode Number (0-8)")
    
    # Save plot as PNG
    plot_filename = os.path.join(output_plot_dir, f"mode_groups_{gp:.2f}h.png")
    plt.savefig(plot_filename, dpi=300)
    plt.close()
    print(f"Saved plot as {plot_filename}")
