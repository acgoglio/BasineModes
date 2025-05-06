import numpy as np
import xarray as xr
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import netCDF4 as ncdf
from matplotlib.colors import LogNorm

mpl.use("Agg")  # For non-interactive backend

# --- Paths ---
indir = "/work/cmcc/ag15419/basin_modes/"
infile = os.path.join(indir, "basin_modes_pow_med.nc")
csvfile = os.path.join(indir, "periods_grouped_greedy_pow.csv")
outfile = os.path.join(indir, "mode_groups_pow.nc")
mesh_mask_file = "/work/cmcc/ag15419/VAA_paper/DATA0/mesh_mask.nc"
bathy_file = "/work/cmcc/ag15419/VAA_paper/DATA0/bathy_meter.nc"
output_plot_dir = os.path.join(indir, "mode_plots_pow")
os.makedirs(output_plot_dir, exist_ok=True)
tolerance = 0.4 # Greedy grouping algorithm (spectrum resolution)

##########################################
# Truncate the colormap to exclude the lightest part (e.g. bottom 20%)
def truncate_colormap(cmap, minval=0.2, maxval=1.0, n=256):
    new_cmap = LinearSegmentedColormap.from_list(
        f"trunc({cmap.name},{minval:.2f},{maxval:.2f})",
        cmap(np.linspace(minval, maxval, n))
    )
    return new_cmap

# --- Load datasets ---
ds = xr.open_dataset(infile)
grouped_df = pd.read_csv(csvfile)
try:
   group_centers = grouped_df["Grouped_Period"].values
except:
   group_centers = grouped_df["Period"].values

# --- Load grid and mask ---
nc_nemogrid = ncdf.Dataset(mesh_mask_file)
nav_lat = nc_nemogrid.variables["nav_lat"][:]
nav_lon = nc_nemogrid.variables["nav_lon"][:]
tmask = nc_nemogrid.variables["tmask"][0, 0, :, :]  # shape (y, x)
field_shape = nav_lat.shape

# Load bathymetry (assumed to be 2D and in meters)
ds_bathy = xr.open_dataset(bathy_file, decode_times=False)
#bathy = ds_bathy["Bathymetry"].values  # Replace "Bathymetry" with correct var name if needed
bathy = ds_bathy["Bathymetry"][0, :, :]

# --- Initialize output fields ---
fields = {f"mode_{gp:.2f}h": np.full(field_shape, np.nan) for gp in group_centers}
pow_fields = {f"pow_{gp:.2f}h": np.full(field_shape, np.nan) for gp in group_centers}

# --- Extract variables ---
modes_vars = [var for var in ds.data_vars if var.startswith("m") and "_T" in var]
pow_vars = {int(var[1]): var for var in ds.data_vars if var.startswith("m") and "_Amp" in var}

# --- Assign values ---
for var in modes_vars:
    mode_num = int(var[1])  # e.g. m0_T -> 0
    period_data = ds[var].values

    # Convert to hours
    if np.issubdtype(period_data.dtype, np.timedelta64):
        period_data = period_data.astype("timedelta64[s]").astype(float) / 3600
    else:
        period_data = period_data.astype(float)
    period_data = np.round(period_data, 2)

    # Get powlitude
    if mode_num not in pow_vars:
        continue
    pow_data = ds[pow_vars[mode_num]].values

    for gp in group_centers:
        match = np.abs(period_data - gp) <= tolerance
        if np.any(match):
            print(f"Mode {mode_num} matches group {gp:.2f}h in {np.sum(match)} points")
        fields[f"mode_{gp:.2f}h"][match] = mode_num
        pow_fields[f"pow_{gp:.2f}h"][match] = np.fmax(pow_fields[f"pow_{gp:.2f}h"][match], pow_data[match])

# --- Save all fields to NetCDF ---
combined_fields = {
    **{name: (("y", "x"), data) for name, data in fields.items()},
    **{name: (("y", "x"), data) for name, data in pow_fields.items()}
}
output_ds = xr.Dataset(
    combined_fields,
    coords={"nav_lat": (("y", "x"), nav_lat), "nav_lon": (("y", "x"), nav_lon)}
)
output_ds.to_netcdf(outfile)
print(f"Saved grouped mode and powlitude maps to: {outfile}")

# --- Plot each grouped mode field ---
for idx_gp,gp in enumerate(group_centers):
    mode_field = fields[f"mode_{gp:.2f}h"]
    pow_field = pow_fields[f"pow_{gp:.2f}h"]

    # Plot 1: mode number
    plt.figure(figsize=(10, 4))
    cmap = mpl.cm.get_cmap("Reds_r")
    cmap.set_bad("white")
    masked_field = np.ma.masked_invalid(mode_field)
    plt.contourf(nav_lon, nav_lat, masked_field, levels=np.arange(0, 9), cmap=cmap, extend="neither")
    plt.colorbar(label="Mode Number")
    plt.contourf(nav_lon, nav_lat, tmask, levels=[-1000, 0.05], colors="gray")
    plt.contour(nav_lon, nav_lat, tmask, levels=[0.05], colors="black", linewidths=0.8)
    #plt.title(f"Modes with Period: {gp:.1f} h ± 0.5 h")
    plt.title(f"Modes with Period: {gp:.2f} h")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.xlim(-6, 36.3)
    plt.ylim(30, 46)
    plt.savefig(os.path.join(output_plot_dir, f"mode_pow_{idx_gp}_{gp:.2f}h.png"), dpi=300, bbox_inches="tight")
    plt.close()

    # Plot 2: presence mask
    presence_mask = (~np.isnan(mode_field)).astype(int)
    plt.figure(figsize=(10, 4))
    cmap_presence = mpl.colors.ListedColormap(["white", "tab:blue"])
    bounds = [0, 0.5, 1.5]
    norm = mpl.colors.BoundaryNorm(bounds, cmap_presence.N)
    plt.contourf(nav_lon, nav_lat, presence_mask, levels=bounds, cmap=cmap_presence, norm=norm)
    plt.contourf(nav_lon, nav_lat, tmask, levels=[-1000, 0.05], colors="gray")
    plt.contour(nav_lon, nav_lat, tmask, levels=[0.05], colors="black", linewidths=0.8)
    #plt.title(f"Presence Map for Period: {gp:.1f} h ± 0.5 h")
    plt.title(f"Presence Map for Period: {gp:.2f} h")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.xlim(-6, 36.3)
    plt.ylim(30, 46)
    # Add legend instead of colorbar
    #legend_elements = [
    #   mpatches.Patch(facecolor="white", edgecolor="black", label="Not present"),
    #   mpatches.Patch(facecolor="tab:blue", edgecolor="black", label="Present")
    #]
    #ax.legend(handles=legend_elements, loc="lower left", frameon=True, fontsize=9)
    plt.savefig(os.path.join(output_plot_dir, f"mode_flag_pow_{idx_gp}_{gp:.2f}h.png"), dpi=300, bbox_inches="tight")
    plt.close()

    # Plot 3: powlitude in % of maximum
    plt.figure(figsize=(10, 4))
    pow_percent = pow_field / np.nanmax(pow_field) * 100
    masked_pow_pct = np.ma.masked_invalid(pow_percent)
    cmap_pow_pct = mpl.cm.get_cmap("gist_stern_r") #("Greens")
    cmap_pow_pct = truncate_colormap(cmap_pow_pct, 0.3, 0.95)
    cmap_pow_pct.set_bad("white")
    im = plt.contourf(nav_lon, nav_lat, masked_pow_pct, levels=np.linspace(0, 100, 41), cmap=cmap_pow_pct)
    plt.colorbar(im, label="Energy (% of max)")

    # Contour at 0%
    plt.contour(nav_lon, nav_lat, pow_percent, levels=[0.005], colors="magenta", linewidths=1.2)

    plt.contourf(nav_lon, nav_lat, tmask, levels=[-1000, 0.05], colors="gray")
    plt.contour(nav_lon, nav_lat, tmask, levels=[0.05], colors="black", linewidths=0.8)

    # Add bathymetry contour lines
    contour_levels = [100, 200, 300]
    plt.contour(nav_lon, nav_lat, bathy, levels=contour_levels, colors="black", linewidths=0.5, linestyles="dashed")
    contour_levels = [400, 500]
    plt.contour(nav_lon, nav_lat, bathy, levels=contour_levels, colors="black", linewidths=0.5, linestyles="dotted")

    #plt.title(f"Energy for Mode with Period: {gp:.1f} h ± 0.5 h")
    plt.title(f"Energy for Mode with Period: {gp:.2f} h")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.xlim(-6, 36.3)
    plt.ylim(30, 46)
    plt.savefig(os.path.join(output_plot_dir, f"mode_powval_{idx_gp}_{gp:.2f}h.png"), dpi=300, bbox_inches="tight")
    plt.close()
