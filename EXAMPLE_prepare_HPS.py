import numpy as np
import obspy, os, gc, re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['image.cmap'] = 'inferno'
plt.rcParams['agg.path.chunksize'] = 10000
from obspy import read, read_inventory

from obspy.imaging.spectrogram import _nearest_pow_2
from scipy import signal
from matplotlib import mlab, transforms

from funda_detection import pitch_detection

## set paths

path_project = '/home/essingd/AWI/knipovich/'
path_waveforms = path_project+'data/'
path_waveforms_test =path_waveforms+'test_folder/'
inv_path = path_project+'CAL/'
    
## set experiment     
# experiment='KNIPA_NEW'
# network = '1B'
# station = 'KNR18'
# component = 'BHZ'


## set preprocessing parameter
pre_filt=(.05,.009,20,24)       # filtercorners for filter during response remove
bandpass_min = .1               # filtercorners for bandpass filter
bandpass_max = 10               # filtercorners for bandpass filter
decimate_factor = 2             # factor to downsample waveform data

## set parameter for spectrum
add_time = 1                    # hours to add before/after day (day to calculate spectrogram) in [hours]; this prevents amplitude drop at beginnig/end of spectrogram
nfft = 2000                     # number of datapoints used in each block for FFT, might be faster with _nearest_pow_2 but due to wanted 80s of win_len in specgram I set it to 2000

## set fundamental frequency detection parameter
nr_downsamp_HPS = 4             # Number to downsample the spectrum for fundamental frequency detection (pitch detection) with HarmonicProductSpectrum-Alogrithm
threshold_HSI = 2               # threshold for ratio of HarmonicStrengthIndex (ratio between amplitude of fundamental_frequency_value and maximum  amplitdue of Interharmonic)
                                # Automated detection and characterization of harmonic tremor in continuous seismic data; Roman D; 2017

trace = read(path_waveforms_test+'1B.KNR08..BHZ.D.2016.325')
trace += read(path_waveforms_test+'1B.KNR08..BHZ.D.2016.326')
trace += read(path_waveforms_test+'1B.KNR08..BHZ.D.2016.327')
trace.merge(method=1) # merges together three days
inv = read_inventory(inv_path+'RESP.1B.KNR08..BHZ')   
print(trace)
## calculate values to trim for 26 hours of data
print(trace)
trace.decimate(factor=decimate_factor) 
trace.remove_response(inventory=inv, pre_filt=pre_filt, output="ACC",water_level=60, plot=False)
trace.filter('bandpass',freqmin=bandpass_min,freqmax=bandpass_max,zerophase=True)
trace.detrend(type='simple') # remove mean
trace.detrend(type='linear') # remove trend
# trace.trim(starttime=trim_start,endtime=trim_end,nearest_sample=True,pad=True, fill_value=0)


    
ff_value_signal, ff_amplitude_signal, ff_value_NOsignal, ff_amplitude_NOsignal, specgram, specgram_time_vec ,specgram_freq_vec,trace_data_one_day,trace_data_one_day_time_line =pitch_detection.pitch_detection(trace,nfft,threshold_HSI,nr_downsamp_HPS)


