import glob
import sys
import xarray as xr
import netCDF4 as nc
import numpy as np
from scipy.signal import periodogram
import heapq
from functools import partial
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy.signal import detrend
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
import shutil
import f_point_powspt
mpl.use('Agg')

########

# Workdir, otufile template and otfil 
work_dir='/work/cmcc/ag15419/basin_modes/'
infile_amppha='/work/cmcc/ag15419/basin_modes/basin_modes_ini.nc'
outfile=work_dir+'basin_modes_pow.nc'

# Infiles
start_date = "20150103"
end_date = "20150201"
all_files = sorted(glob.glob("/work/cmcc/ag15419/exp/fix_mfseas9_longrun_hmslp_2NT_AB/EXP00_BF/20*/model/medfs-eas9_1h_20*_2D_grid_T.nc"))

# Model time step in seconds
dt = 3600

###################A
# Build the outfile:
shutil.copy(infile_amppha,outfile)

# Read infiles
# Select the period
infile = []
for f in all_files:
    parts = f.split("/")
    file_date = parts[7] # 7 6  
    if start_date <= file_date <= end_date: 
            infile.append(f)

# Read lat and lon in the first file
nav_lat = None
nav_lon = None
first=0
for nc2open in infile:
  if first==0:
    model = nc.Dataset(nc2open, 'r')
    if nav_lat is None:  
        nav_lat = model.variables['nav_lat'][:]
        nav_lon = model.variables['nav_lon'][:]
        first=1
    model.close()

# Read ssh xarray
ds = xr.open_mfdataset(infile, combine='by_coords', parallel=True)
ssh_ts_all = ds['sossheig']

########################
# Compute and write values in the netCDF file
modes_outfile = nc.Dataset(outfile, 'a')

# Call the function for each point in the Med
for lon_idx in range (300,len(nav_lon)): #(0,len(nav_lon)):
    for lat_idx in range (0,len(nav_lat)): # (0,len(nav_lat)):
        ssh_ts_point = ssh_ts_all[:, lat_idx, lon_idx].values
        pow_peak_periods_main, pow_peak_amplitudes_main = f_point_powspt.pow_main_modes(lat_idx, lon_idx, ssh_ts_point, dt)

        for i in range(8):
            try:
                #print(f'Mode {i} T={amp_peak_periods_main[i]:.2f} h, Amp={amp_peak_amplitudes_main[i]:.4f} m')
                modes_outfile.variables[f'm{i}_Amp'][lat_idx, lon_idx] = pow_peak_amplitudes_main[i]
                modes_outfile.variables[f'm{i}_T'][lat_idx, lon_idx] = pow_peak_periods_main[i]
            except:
                modes_outfile.variables[f'm{i}_Amp'][lat_idx, lon_idx] = np.nan
                modes_outfile.variables[f'm{i}_T'][lat_idx, lon_idx] = np.nan
modes_outfile.close()
