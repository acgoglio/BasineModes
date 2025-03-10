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
mpl.use('Agg')

########
# Inputs and outputs

start_date = "20150201"
end_date = "20150601"

#all_files = sorted(glob.glob("/work/cmcc/ag15419/exp/fix_mfseas9_longrun_barotropic_final22/EXP00/2*/model/medfs-eas9_1h_2*_2D_grid_T.nc"))
#all_files = sorted(glob.glob("/work/cmcc/ag15419/exp/fix_mfseas9_longrun_surge_2/EXP00/2016*/model/medfs-eas9_1h_2016*_2D_grid_T.nc"))
#all_files = sorted(glob.glob("/work/cmcc/med-dev/exp/EAS9_assw_nt/202*/model/medfs-eas9_1h_202*_2D_grid_T.nc"))
all_files = sorted(glob.glob("/work/cmcc/ag15419/exp/fix_mfseas9_longrun_surge_2NT/EXP00_2s/20*/model/medfs-eas9_1ts_2015*_2D_grid_T.nc"))

# Exp tag
exp='4m_NT'

# Lat and lon indexes
lat_idx = str(sys.argv[1]) #138 #358 #360
lon_idx = str(sys.argv[2]) #331 #744 #746

# Outfile
outfile = str(sys.argv[3])

# Model time step in seconds
dt = 3*60

# Number of modes to analyze
n_modes = 4  

# Minimum peak amplitude ; min,max width ; min distance between peaks to detect peaks in Amp plots (meters,hours, points respectively)
amp_peak_height=0.0001
amp_peak_width=(0, 40)
amp_peak_distance=3

# Flag and threshold [h] for filtering the spectrum the threshold is also used as plot minimum 
flag_filter='true'
th_filter=72

# Flag for spectrum detrending:
flag_detrend='true'

# Flag for Gaussian smoothing of the spectrum
flag_smooth='true'
sigma=4
#def moving_average(data, window_size):
#    return np.convolve(data, np.ones(window_size) / window_size, mode='same')
#window_size=11

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
    #print('Processing:', nc2open)
    model = nc.Dataset(nc2open, 'r')
    ssh_ts = np.array(model.variables['sossheig'][:, lat_idx, lon_idx])
    ssh_ts_all = np.concatenate((ssh_ts_all, ssh_ts))

    if not grid_info:
        lats = model.variables['nav_lat'][lat_idx, lon_idx]
        lons = model.variables['nav_lon'][lat_idx, lon_idx]
        grid_info = True

    #print(f"Total number of points in the SSH time series: {len(ssh_ts_all)}")
    model.close()

# Convert SSH time series to NumPy array and remove NaNs
time_series_point = np.array(ssh_ts_all)
#print ('Check input time series:',time_series_point.shape,time_series_point)
valid_indices = np.logical_not(np.isnan(time_series_point))
time_series_clean = time_series_point[valid_indices]
#print ('Check clean time series:',time_series_clean.shape,time_series_clean)

#### SPECTRUM ANALYSIS

# Compute FFT
spt_len = len(time_series_clean)
#print('Time series values:', spt_len)
fft = np.fft.fft(time_series_clean)
freq = np.fft.fftfreq(spt_len, d=dt)

# Select only positive frequencies (first half of FFT, excluding Nyquist frequency)
half_spt_len = spt_len // 2
freq_positive = freq[:half_spt_len]  # Only positive frequencies
fft_positive = fft[:half_spt_len]  # Only corresponding FFT values

# Ensure we are not including zero frequency
freq_positive = freq_positive[freq_positive > 0]
fft_positive = fft_positive[:len(freq_positive)]  # Match the length of fft_positive with freq_positive

# Compute Power Spectrum Density (normalized)
spt = (np.abs(fft_positive) ** 2) / spt_len**2

# Compute Periods in hours
periods = 1 / freq_positive / 3600  # Convert periods to hours
# Check that periods span the expected range
#print(f"Periods (hours) range from {periods[0]:.2f} h to {periods[-1]:.2f} h.")

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
   spt_smooth = gaussian_filter1d(spt, sigma=sigma)
   #spt_smooth = moving_average(spt, window_size)
   amp_smooth = gaussian_filter1d(amplitudes, sigma=sigma)
else:
   spt_smooth = spt
   amp_smooth = amplitudes

# Detrend
if flag_detrend == 'true':
   spt_det = detrend(spt_smooth)
   #detrend_window_size = int(spt_len/2)
   #for i in range(0, len(spt_smooth), detrend_window_size):
   #    spt_det[i:i+window_size] = detrend(spt_smooth[i:i+window_size])
else:
   spt_det = spt

# Found peaks in spt and in amplitude:
peaks, _ = find_peaks(spt_smooth,height=0.000001,distance=5)
peak_frequencies = freq_positive[peaks]
#peak_amplitudes = spt_smoothed[peaks]

amp_peaks, _ = find_peaks(amp_smooth,prominence=amp_peak_height,width=amp_peak_width,distance=amp_peak_distance)
amp_peak_frequencies = freq_positive[amp_peaks]
amp_peak_amplitudes = amp_smooth[amp_peaks]

# Order based on amplitudes
sorted_indices_peak_amp = np.argsort(amp_peak_amplitudes)[::-1]
amp_peak_amplitudes_sorted = amp_peak_amplitudes[sorted_indices_peak_amp]
amp_peak_frequencies_sorted = amp_peak_frequencies[sorted_indices_peak_amp]

# Time array for SSH plot (convert to hours)
ssh_time = np.arange(0, spt_len * dt, dt) / 3600  

# Select the main spectral peaks based on spectral density
top_indices_spt_det = np.argsort(spt_det)[-n_modes:]  # Get indices of top n_modes power values
sorted_indices_spt_det = np.argsort(spt_det[top_indices_spt_det])[::-1]  # Sort by descending power

# Extract the corresponding frequencies, periods, and amplitudes for density-based selection
top_freq_positive_spt_det = freq_positive[top_indices_spt_det][sorted_indices_spt_det]
top_periods_spt_det = periods[top_indices_spt_det][sorted_indices_spt_det]
top_amplitudes_spt_det = amplitudes[top_indices_spt_det][sorted_indices_spt_det]

# Print the main spectral modes based on spectral density
#print("Main Spectral Modes (based on power spectrum density):")
#for i in range(n_modes):
#    print(f"Mode {i+1}: Period = {top_periods_spt_det[i]:.2f} h, Freq = {top_freq_positive_spt_det[i]:.6f} Hz, Amp = {top_amplitudes_spt_det[i]:.3f} m")

# Now select the main modes based on amplitude
n_valid = min(len(amplitudes), len(freq_positive))  # Ensure valid index range
top_indices_amp = np.argpartition(amplitudes[:n_valid], -n_modes)[-n_modes:]  # Select indices of top amplitudes
sorted_indices_amp = np.argsort(amplitudes[top_indices_amp])[::-1]  # Sort by descending amplitude

# Extract the corresponding frequencies, periods, and amplitudes for amplitude-based selection
top_freq_positive_amp = freq_positive[top_indices_amp][sorted_indices_amp]
top_periods_amp = periods[top_indices_amp][sorted_indices_amp]
top_amplitudes_amp = amplitudes[top_indices_amp][sorted_indices_amp]

# Print the main spectral modes based on amplitude
#print("Main Spectral Modes (based on amplitude):")
#for i in range(n_modes):
#    print(f"Mode {i+1}: Period = {top_periods_amp[i]:.2f} h, Freq = {top_freq_positive_amp[i]:.6f} Hz, Amp = {top_amplitudes_amp[i]:.3f} m")

#######################
# PLOT SSH
plt.figure(figsize=(10, 6))
plt.title(f'SSH at lat={lats} lon={lons}')
plt.plot(ssh_time, time_series_clean, '-', label=f'SSH at lat={lats} lon={lons}')
plt.xlabel('Time (h)')
plt.ylabel('SSH (m)')
plt.grid()
plt.legend()
plt.savefig(f'ssh_{lat_idx}_{lon_idx}_{exp}.png')

# PLOT POWER SPECTRUM
#plt.figure(figsize=(15, 9))
#plt.title(f'Power Spectrum at lat={lats} lon={lons}')
#plt.loglog(periods, spt, marker='o', linestyle='-', label='Power Spectrum')
#if flag_smooth == 'true':
#   plt.loglog(periods, spt_smooth, marker='o', linestyle='-', label='Smoothed Power Spectrum')
#plt.xlabel('Period (h)')
#plt.ylabel('Power Spectrum')
#plt.axvline(24, color='black', linestyle='-')
#plt.axvline(12, color='black', linestyle='-')
#plt.axvline(6, color='black', linestyle='-')
#
## Mark the main modes based on spectral density
#for i in range(n_modes):
#    plt.axvline(top_periods_spt_det[i], color='red', linestyle='--', label=f'SPT Mode {i+1} (T = {top_periods_spt_det[i]:.2f} h, Amp = {top_amplitudes_spt_det[i]:.3f} m)')
#
## Mark the main modes based on amplitude
#for i in range(n_modes):
#    plt.axvline(top_periods_amp[i], color='blue', linestyle='--', label=f'Amp Mode {i+1} (T = {top_periods_amp[i]:.2f} h, Amp = {top_amplitudes_amp[i]:.3f} m)')
#
## Mark the main modes based on peak finder
#for i in range(0,len(peak_frequencies)):
#    plt.axvline(1/peak_frequencies[i]/3600, color='green',linestyle='--')
#
#plt.xlim(th_filter,0.5)
#plt.grid()
#plt.legend()
#plt.savefig(f'spt_{lat_idx}_{lon_idx}_{exp}.png')
#
## Non-log spt plot
#plt.figure(figsize=(15, 9))
#plt.title(f'Power Spectrum at lat={lats} lon={lons}')
#plt.loglog(periods, spt, marker='o', linestyle='-', label='Power Spectrum')
#if flag_smooth == 'true':
#   plt.loglog(periods, spt_smooth, marker='o', linestyle='-', label='Smoothed Power Spectrum')
#plt.xlabel('Period (h)')
#plt.ylabel('Power Spectrum')
#plt.axvline(24, color='black', linestyle='-')
#plt.axvline(12, color='black', linestyle='-')
#plt.axvline(6, color='black', linestyle='-')
#
## Mark the main modes based on spectral density
#for i in range(n_modes):
#    plt.axvline(top_periods_spt_det[i], color='red', linestyle='--', label=f'SPT Mode {i+1} (T = {top_periods_spt_det[i]:.2f} h, Amp = {top_amplitudes_spt_det[i]:.3f} m)')
#
## Mark the main modes based on amplitude
#for i in range(n_modes):
#    plt.axvline(top_periods_amp[i], color='blue', linestyle='--', label=f'Amp Mode {i+1} (T = {top_periods_amp[i]:.2f} h, Amp = {top_amplitudes_amp[i]:.3f} m)')
#
#plt.xlim(th_filter,0.5)
#plt.yscale('linear')
#plt.grid()
#plt.legend()
#plt.savefig(f'spt_nolog_{lat_idx}_{lon_idx}_{exp}.png')
#
## Peaks plot
#plt.figure(figsize=(15, 9))
#plt.title(f'Detrended Power Spectrum at lat={lats} lon={lons}')
#plt.loglog(periods, spt_det, marker='o', linestyle='-', label='Detrended Power Spectrum')
##if flag_smooth == 'true':
##   plt.loglog(periods, spt_smooth, marker='o', linestyle='-', label='Smoothed Power Spectrum')
#plt.xlabel('Period (h)')
#plt.ylabel('Power Spectrum')
#plt.axvline(24, color='black', linestyle='-')
#plt.axvline(12, color='black', linestyle='-')
#plt.axvline(6, color='black', linestyle='-')
#
## Mark the main modes based on spectral density
#for i in range(n_modes):
#    plt.axvline(top_periods_spt_det[i], color='red', linestyle='--', label=f'SPT Mode {i+1} (T = {top_periods_spt_det[i]:.2f} h, Amp = {top_amplitudes_spt_det[i]:.3f} m)')
#
## Mark the main modes based on amplitude
#for i in range(n_modes):
#    plt.axvline(top_periods_amp[i], color='blue', linestyle='--', label=f'Amp Mode {i+1} (T = {top_periods_amp[i]:.2f} h, Amp = {top_amplitudes_amp[i]:.3f} m)')
#
#plt.xlim(th_filter,0.5)
#plt.grid()
#plt.legend()
#plt.savefig(f'spt_det_{lat_idx}_{lon_idx}_{exp}.png')

# Amplitude plots
plt.figure(figsize=(15, 9))
plt.title(f'Modes amplitudes at lat={lats} lon={lons}')
plt.loglog(periods, amplitudes, marker='o', linestyle='-', label='Modes Amplitudes')
plt.loglog(periods, amp_smooth, marker='o', linestyle='-', label='Smoothed Modes Amplitudes')
plt.xlabel('Period (h)')
plt.ylabel('Mode Amplitude (m)')
plt.axvline(24, color='black', linestyle='-')
plt.axvline(12, color='black', linestyle='-')
plt.axvline(6, color='black', linestyle='-')

# Mark the main modes based on peak finder
for i in range(0,len(amp_peak_frequencies)):
    plt.axvline(1/amp_peak_frequencies[i]/3600, color='green',linestyle='--',label=f'Mode {i} (T={1/amp_peak_frequencies_sorted[i]/3600:.2f} h, Amp={amp_peak_amplitudes_sorted[i]:.4f} m)')
plt.xlim(th_filter-1,0.5)
plt.grid()
plt.legend()
plt.savefig(f'amp_{lat_idx}_{lon_idx}_{exp}.png')

# Amp no log plot
#plt.figure(figsize=(15, 9))
#plt.title(f'Modes amplitudes at lat={lats} lon={lons}')
#plt.loglog(periods, amplitudes, marker='o', linestyle='-', label='Power Spectrum')
#plt.xlabel('Period (h)')
#plt.ylabel('Mode Amplitude (m)')
#plt.axvline(24, color='black', linestyle='-')
#plt.axvline(12, color='black', linestyle='-')
#plt.axvline(6, color='black', linestyle='-')

#plt.xlim(th_filter,0.5)
#plt.yscale('linear')
#plt.grid()
#plt.legend()
#plt.savefig(f'amp_nolog_{lat_idx}_{lon_idx}_{exp}.png')


########################
# Write values in the netCDF file
out_lat_idx=0
out_lon_idx=0
modes_outfile = nc.Dataset(outfile, 'a')
for i in range(0,8):
       var_amp = modes_outfile.variables['m'+str(i)+'_Amp']
       var_T = modes_outfile.variables['m'+str(i)+'_T']
       try:
          print (f'Mode {i} T={1/amp_peak_frequencies_sorted[i]/3600:.2f} h, Amp={amp_peak_amplitudes_sorted[i]:.4f} m')
          var_amp[out_lat_idx,out_lon_idx] = amp_peak_amplitudes_sorted[i]
          var_T[out_lat_idx,out_lon_idx] = 1/amp_peak_frequencies_sorted[i]/3600
       except:
          var_amp[out_lat_idx,out_lon_idx] = np.nan
          var_T[out_lat_idx,out_lon_idx] = np.nan

modes_outfile.close()

