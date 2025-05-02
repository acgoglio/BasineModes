import numpy as np
import xarray as xr
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import netCDF4 as ncdf

mpl.use("Agg")  # For non-interactive backend

# --- Paths ---
indir = "/work/cmcc/ag15419/basin_modes/"
infile = os.path.join(indir, "basin_modes_pow_med.nc")
csvfile = os.path.join(indir, "periods_grouped30min_pow.csv")
outfile = os.path.join(indir, "mode_groups_pow.nc")
mesh_mask_file = "/work/cmcc/ag15419/VAA_paper/DATA0/mesh_mask.nc"
output_plot_dir = os.path.join(indir, "mode_plots_pow")
os.makedirs(output_plot_dir, exist_ok=True)

# --- Load datasets ---
ds = xr.open_dataset(infile)
grouped_df = pd.read_csv(csvfile)
group_centers = grouped_df["Grouped_Period"].values

# --- Load grid and mask ---
nc_nemogrid = ncdf.Dataset(mesh_mask_file)
nav_lat = nc_nemogrid.variables["nav_lat"][:]
nav_lon = nc_nemogrid.variables["nav_lon"][:]
tmask = nc_nemogrid.variables["tmask"][0, 0, :, :]  # shape (y, x)
field_shape = nav_lat.shape

# --- Initialize output fields with np.nan ---
fields = {f"mode_{gp:.2f}h": np.full(field_shape, np.nan) for gp in group_centers}

# --- Extract modes and assign values ---
modes_vars = [var for var in ds.data_vars if var.startswith("m") and "_T" in var]

for var in modes_vars:
    mode_num = int(var[1])  # Extract the integer from variable name like m0_T
    mode_period_raw = ds[var].values

    if np.issubdtype(mode_period_raw.dtype, np.timedelta64):
        mode_period = mode_period_raw.astype("timedelta64[s]").astype(float) / 3600
    else:
        mode_period = mode_period_raw.astype(float)

    mode_period = np.round(mode_period, 2)

    for gp in group_centers:
        match = np.abs(mode_period - gp) <= 0.5
        if np.any(match):
            print(f"Mode {mode_num} matches group {gp:.2f}h in {np.sum(match)} points")
        current_field = fields[f"mode_{gp:.2f}h"]
        current_field[match] = mode_num

# --- Save to NetCDF ---
output_ds = xr.Dataset(
    {name: (("y", "x"), data) for name, data in fields.items()},
    coords={"nav_lat": (("y", "x"), nav_lat), "nav_lon": (("y", "x"), nav_lon)}
)
output_ds.to_netcdf(outfile)
print(f"Saved grouped mode maps to: {outfile}")

# --- Plot each grouped mode field ---
for gp in group_centers:
    mode_field = fields[f"mode_{gp:.2f}h"]

    #################################
    # Plot the mod num per resonance period
    plt.figure(figsize=(10, 4))

    cmap = mpl.cm.get_cmap("Reds_r") #.copy()
    cmap.set_bad("white")

    # Mask nan for plotting
    masked_field = np.ma.masked_invalid(mode_field)

    plt.contourf(nav_lon, nav_lat, masked_field, levels=np.arange(0, 9), cmap=cmap, extend="neither")
    plt.colorbar(label="Mode Number")

    # Add land-sea contour
    plt.contourf(nav_lon, nav_lat, tmask, levels=[-1000,0.05], colors="gray")
    plt.contour(nav_lon, nav_lat, tmask, levels=[0.5], colors="black", linewidths=0.8)

    plt.title(f"Modes with Period: {gp:.1f} h +- 0.5 h")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.xlim(-6, 36.3)
    plt.ylim(30, 46)

    plot_file = os.path.join(output_plot_dir, f"mode_pow_{gp:.1f}h.png")
    plt.savefig(plot_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved plot as {plot_file}")

    #################################
    # Plot the flag per resonance period
    # --- Second plot: highlight area with data (non-NaN) in red ---
    # --- Second plot: highlight area with data (non-NaN) in red ---
    presence_mask = (~np.isnan(mode_field)).astype(int)  # 1 where data is present, 0 elsewhere

    plt.figure(figsize=(10, 4))

    cmap_presence = mpl.colors.ListedColormap(["white", "tab:blue"])
    bounds = [0, 0.5, 1.5]
    norm = mpl.colors.BoundaryNorm(bounds, cmap_presence.N)

    plt.contourf(nav_lon, nav_lat, presence_mask, levels=bounds, cmap=cmap_presence, norm=norm)

    # Add land-sea contour
    plt.contourf(nav_lon, nav_lat, tmask, levels=[-1000,0.05], colors="gray")
    plt.contour(nav_lon, nav_lat, tmask, levels=[0.05], colors="black", linewidths=0.8)

    plt.title(f"Mode with Period {gp:.1f} h +- 0.5 h")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.xlim(-6, 36.3)
    plt.ylim(30, 46)

    presence_plot_file = os.path.join(output_plot_dir, f"mode_flag_pow_{gp:.1f}h.png")
    plt.savefig(presence_plot_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved presence map as {presence_plot_file}")


