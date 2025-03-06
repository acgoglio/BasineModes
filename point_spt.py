import glob
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
mpl.use('Agg')

########
# Inputs and outputs

# Lat and lon indexes
lat_idx = argv[0]
lon_idx = argv[1]

# Outfile
outfile = argv[2]

start_date = "20150201"
end_date = "20150601"

all_files = sorted(glob.glob("/work/cmcc/ag15419/exp/fix_mfseas9_longrun_surge_2NT/EXP00_2s/20*/model/medfs-eas9_1ts_2015*_2D_grid_T.nc"))

# Model time step in seconds
dt = 3*60

# Number of modes to analyze
n_modes = 4  

# Minimum peak amplitude ; min,max width ; min distance between peaks to detect peaks in Amp plots (meters,hours, points respectively)
amp_peak_height=0.0005
amp_peak_width=(0, 100)
amp_peak_distance=4

# Flag and threshold [h] for filtering the spectrum the threshold is also used as plot minimum 
flag_filter='true'
th_filter=72

# Flag for spectrum detrending:
flag_detrend='true'

# Flag for Gaussian smoothing of the spectrum
flag_smooth='true'
sigma=4

###################
# Select the period
infile = []
for f in all_files:
    parts = f.split("/")
    file_date = parts[7] # 7 6  
    if start_date <= file_date <= end_date: 
            infile.append(f)

# Initialize SSH time series
ssh_ts_all = []

# Read data from NetCDF files
grid_info = False
for nc2open in infile:
    model = nc.Dataset(nc2open, 'r')
    ssh_ts = np.array(model.variables['sossheig'][:, lat_idx, lon_idx])
    ssh_ts_all = np.concatenate((ssh_ts_all, ssh_ts))

    if not grid_info:
        lats = model.variables['nav_lat'][lat_idx, lon_idx]
        lons = model.variables['nav_lon'][lat_idx, lon_idx]
        grid_info = True

    model.close()

# Convert SSH time series to NumPy array and remove NaNs
time_series_point = np.array(ssh_ts_all)
valid_indices = np.logical_not(np.isnan(time_series_point))
time_series_clean = time_series_point[valid_indices]

#### SPECTRUM ANALYSIS

# Compute FFT
spt_len = len(time_series_clean)
fft = np.fft.fft(time_series_clean)
freq = np.fft.fftfreq(spt_len, d=dt)

# Select only positive frequencies (first half of FFT, excluding Nyquist frequency)
half_spt_len = spt_len // 2
freq_positive = freq[:half_spt_len]  # Only positive frequencies
fft_positive = fft[:half_spt_len]  # Only corresponding FFT values

# Ensure we are not including zero frequency
freq_positive = freq_positive[freq_positive > 0]
fft_positive = fft_positive[:len(freq_positive)]  # Match the length of fft_positive with freq_positive

# Compute Periods in hours
periods = 1 / freq_positive / 3600  # Convert periods to hours

# Compute Amplitudes (correct scaling)
amplitudes = (2 / spt_len) * np.abs(fft_positive)

if flag_filter=='true':
   # Apply high-pass filter: Set frequencies below the threshold to zero
   high_pass_threshold = 1 / (th_filter * 3600)  # Corresponding to 48 hours in Hz
   fft_positive[freq_positive < high_pass_threshold] = 0  # Filter out frequencies below the threshold

   # Recompute the power spectrum and amplitudes after filtering
   spt = (np.abs(fft_positive) ** 2) / spt_len**2
   amplitudes = (2 / spt_len) * np.abs(fft_positive)

# Smooth
if flag_smooth == 'true':
   amp_smooth = gaussian_filter1d(amplitudes, sigma=sigma)
else:
   amp_smooth = amplitudes

# Found peaks in spt and in amplitude:

amp_peaks, _ = find_peaks(amp_smooth,prominence=amp_peak_height,width=amp_peak_width,distance=amp_peak_distance)
amp_peak_frequencies = freq_positive[amp_peaks]
amp_peak_amplitudes = amp_smooth[amp_peaks]

# Time array for SSH plot (convert to hours)
ssh_time = np.arange(0, spt_len * dt, dt) / 3600  

########################
# Write values in the netCDF file
modes_outfile = nc.Dataset(outfile, 'r+')
for i in range(0,7):
       var_amp = modes_outfile.variables['m'+i+'_Amp']
       var_T = modes_outfile.variables['m'+i+'_T']
       try:
          #print (f'Mode {i} T={1/amp_peak_frequencies[i]/3600:.2f} h, Amp={amp_peak_amplitudes[i]:.3f} m')
          var_amp[lat_idx, lon_idx] = amp_peak_amplitudes[i]
          var_T[lat_idx, lon_idx] = 1/amp_peak_frequencies[i]/3600
       except:
          var_amp[lat_idx, lon_idx] = np.nan
          var_T[lat_idx, lon_idx] = np.nan

modes_outfile.close()
