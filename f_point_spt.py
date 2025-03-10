
def# -------------------------------------------
# Extract the main 8 modes at each grid point
# -------------------------------------------

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

def amp_main_modes(lat_idx,lon_idx,ssh_ts_all,dt):

   ########
   
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
   
   ###################
   
   # Convert SSH time series to NumPy array and remove NaNs
   time_series_point = np.array(ssh_ts_all)
   valid_indices = np.logical_not(np.isnan(time_series_point))
   time_series_clean = time_series_point[valid_indices]
   
    if len(time_series_clean) == 0:
        return np.full(8, np.nan), np.full(8, np.nan)

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
   
   # Found peaks in amplitude:
   amp_peaks, _ = find_peaks(amp_smooth,prominence=amp_peak_height,width=amp_peak_width,distance=amp_peak_distance)
   amp_peak_frequencies = freq_positive[amp_peaks]
   amp_peak_amplitudes = amp_smooth[amp_peaks]
   
   # Order based on amplitudes
   sorted_indices = np.argsort(amp_peak_amplitudes)[::-1]
   amp_peak_amplitudes_sorted = amp_peak_amplitudes[sorted_indices]
   amp_peak_frequencies_sorted = amp_peak_frequencies[sorted_indices]
   
   amp_peak_periods_main=[]
   amp_peak_amplitudes_main=[]
   for i in range(0,8):
     try:
       amp_peak_periods_main.append(1/amp_peak_frequencies_sorted[i]/3600)
       amp_peak_amplitudes_main.append(amp_peak_amplitudes_sorted[i])
     except:
       amp_peak_periods_main.append(np.nan)
       amp_peak_amplitudes_main.append(np.nan)
   
   return np.array(amp_peak_periods_main),np.array(amp_peak_amplitudes_main)

