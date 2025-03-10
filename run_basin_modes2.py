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
import f_point_spt
mpl.use('Agg')

########

# Workdir, otufile template and otfil 
work_dir='/work/cmcc/ag15419/basin_modes/'
infile_amppha='/work/cmcc/ag15419/basin_modes/basin_modes_ini.nc'
outfile=work_dir+'basin_modes.nc'

# Infiles
start_date = "20150201"
end_date = "20150601"
all_files = sorted(glob.glob("/work/cmcc/ag15419/exp/fix_mfseas9_longrun_surge_2NT/EXP00_2s/20*/model/medfs-eas9_1ts_2015*_2D_grid_T.nc"))
# Model time step in seconds
dt = 3*60

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

# Initialize SSH time series
#ssh_ts_list = []  
nav_lat = None
nav_lon = None

first=0
for nc2open in infile:
    model = nc.Dataset(nc2open, 'r')
    #ssh_ts_list.append(model.variables['sossheig'][:])

    if nav_lat is None:  # First file: read lat/lon
        nav_lat = model.variables['nav_lat'][:]
        nav_lon = model.variables['nav_lon'][:]

    model.close()

# Cat all files
#ssh_ts_all = np.concatenate(ssh_ts_list, axis=0)

# Read ssh xarray
ds = xr.open_mfdataset(infile, combine='by_coords', parallel=True)
ssh_ts_all = ds['sossheig']

########################
# Compute and write values in the netCDF file
modes_outfile = nc.Dataset(outfile, 'a')
#for i in range(0,8):
#       var_amp = modes_outfile.variables['m'+str(i)+'_Amp']
#       var_T = modes_outfile.variables['m'+str(i)+'_T']

# Call the function for each point
for lon_idx in range (285,len(nav_lon)): #(0,len(nav_lon)):
    for lat_idx in range (0,len(nav_lat)): # (0,len(nav_lat)):
        ssh_ts_point = ssh_ts_all[:, lat_idx, lon_idx].values
        amp_peak_periods_main, amp_peak_amplitudes_main = f_point_spt.amp_main_modes(lat_idx, lon_idx, ssh_ts_point, dt)
        #amp_peak_periods_main,amp_peak_amplitudes_main=f_point_spt.amp_main_modes(lat_idx,lon_idx,ssh_ts_all[:, lat_idx, lon_idx],dt)

        for i in range(8):  # Ci sono 8 modi
            try:
                print(f'Mode {i} T={amp_peak_periods_main[i]:.2f} h, Amp={amp_peak_amplitudes_main[i]:.4f} m')
                modes_outfile.variables[f'm{i}_Amp'][lat_idx, lon_idx] = amp_peak_amplitudes_main[i]
                modes_outfile.variables[f'm{i}_T'][lat_idx, lon_idx] = amp_peak_periods_main[i]
            except:
                modes_outfile.variables[f'm{i}_Amp'][lat_idx, lon_idx] = np.nan
                modes_outfile.variables[f'm{i}_T'][lat_idx, lon_idx] = np.nan
              #try:
              #   print (f'Mode {i} T={amp_peak_periods_main[i]} h, Amp={amp_peak_amplitudes_main[i]} m')
              #   var_amp[lat_idx,lon_idx] = amp_peak_amplitudes_main[i]
              #   var_T[lat_idx,lon_idx] = amp_peak_periods_main[i]
              #except:
              #   var_amp[out_lat_idx,out_lon_idx] = np.nan
              #   var_T[out_lat_idx,out_lon_idx] = np.nan

modes_outfile.close()
