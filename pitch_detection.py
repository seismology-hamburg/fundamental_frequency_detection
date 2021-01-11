
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2020-12-15T21-00-00

@author: dessing
"""


import obspy
import numpy as np
import numpy.ma as ma
from matplotlib import mlab, transforms




def hps_algo(input_data, freq_vector, nr_downsamp,min_freq=0, max_freq=0.5):
    '''HarmonicProductSpectrum Algorithum to obtain fundamental frequency  
    
    http://musicweb.ucsd.edu/~trsmyth/analysis/Harmonic_Product_Spectrum.html 
    
    INPUT
    
    input_data                          ndarray      spectrogram including harmonic signal
    freq_vector                         ndarray      frequency vector of spectrogram
    min_freq                            float        minimal frequency to mask frequency band of microseismic noise, optional
    max_freq                            float        maximal frequency to mask frequency band of microseismic noise, optional
    nr_dwonsamp                         integer      value to downsamp 
    
    OUTPUT
    ff_value                            ndarray     value of detected fundamental frequency (FF)
    ff_amplitude                        ndarray     amplitude values for FF
    iharmonic_ival_max_amp              ndarray     maximum value found in the interharmonic (IH) part
    
    
       FF1------IH-----FF2------IH----FF3
        
        |               |              |       
        |               |              |
        |               |              |
        |_______________|______________|
        -----------------------------------> 
                Frequency
    
    '''
    
    
    # mask array in given frequency band
    freq = freq_vector
    
    lower_minval_freq = min_freq
    lower_maxval_freq = max_freq
    lower_freq_values = freq[np.where(np.logical_and(freq>lower_minval_freq, freq<lower_maxval_freq))]

    lower_first_indi_to_use = np.where(freq == lower_freq_values[0])[0]
    lower_last_indi_to_use = np.where(freq == lower_freq_values[-1])[0]
    input_data_masked = ma.array(input_data)
    input_data_masked[int(lower_first_indi_to_use):int(lower_last_indi_to_use),:] = 0
    harmonic_product_spectrum = np.zeros_like(input_data)
    
    # 1D arrays
    ff_value = np.zeros((np.shape(input_data))[1])
    ff_indi = np.zeros_like(ff_value)
    ff_amplitude = np.zeros_like(ff_value)
    iharmonic_ival_max_amp = np.zeros_like(ff_value)

    
    spectrum_downsamp =  np.zeros(((np.shape(input_data_masked))[1],nr_downsamp, len(input_data_masked)))
    
    
    # for loop over all columns
    for col in range((np.shape(input_data_masked))[1]):
        spectrum_signal = input_data_masked[:,col]
        spectrum_downsamp[col,0,:] = spectrum_signal
    
        # for loop to downsample every column in spectrogram
        for i in range(nr_downsamp-1):
            nr_downsamp_array = np.linspace(2,nr_downsamp, nr_downsamp-1)
            spectrum_downsamp[col,i+1,:][0:len(spectrum_signal[::int(nr_downsamp_array[i])])] = spectrum_signal[::int(nr_downsamp_array[i])]

        if nr_downsamp==4:
            harmonic_product_spectrum[:,col] = spectrum_downsamp[col,0,:] *  spectrum_downsamp[col,1,:] * spectrum_downsamp[col,2,:] * spectrum_downsamp[col,3,:]

        if nr_downsamp==3:
            harmonic_product_spectrum[:,col] = spectrum_downsamp[col,0,:] *  spectrum_downsamp[col,1,:] * spectrum_downsamp[col,2,:]

        if nr_downsamp==2:
            harmonic_product_spectrum[:,col] = spectrum_downsamp[col,0,:] *  spectrum_downsamp[col,1,:]
            
        ff_value[col] = freq[np.argmax(harmonic_product_spectrum[:,col])]
        ff_indi[col] = np.argmax(harmonic_product_spectrum[:,col])
        ff_amplitude[col] = spectrum_signal[int(ff_indi[col])]
        iharmonic_ival = spectrum_signal[int(ff_indi[col]*1.25):int(ff_indi[col]*1.75)]
        iharmonic_ival_max_amp[col] = np.max(iharmonic_ival)
       
            
    return ff_value, ff_amplitude, iharmonic_ival_max_amp 





def pitch_detection(trace,nfft,threshold_HSI,nr_downsamp_HPS,add_time=1):
    '''Pitch deteciton for python after the idea of Roman, D. C. (2017), Automated detection and characterization of harmonic tremor in continuous seismic data, Geophys. Res. Lett., 44, 6065–6073 doi:10.1002/2017GL073715
    

    
    INPUT
    
    trace                               obspy.trace  time domain signal of three days of signal in a row, the day in the center will be checked for harmonic tremor
    nfft                                integer      frequency vector of spectrogram
    threshold_HSI                       integer      threshold for Harmonic Strength Index (HSI) which is the ration between the amplitude of the fundamental frequency and the interharmonic range
    nr_downsamp_HPS                     integer      value to downsample the spectrum for fundamental frequency detection
    add_time                            integer      add time (in hours) before and after day of interest to hinder sharp drop at edge values 
    
    OUTPUT 
    ff_value_signal                     ndarray     fundamental frequency values of detected harmonic tremor 
    ff_amplitude_signal                 ndarray     amplitude values of detected harmonic tremor
    ff_value_NOsignal                   ndarray     fundamental frequency values rejected by HSI
    ff_amplitude_NOsignal               ndarray     amplitude values rejected by HSI
    specgram                            ndarray     spectrogram calculated 
    specgram_time_vec                   ndarray     time vector of spectrogram
    specgram_freq_vec                   ndarray     frequency vector of spectrogram
    trace_data_one_day                  ndarray     time domain signal of day of interest
    trace_data_one_day_time_line        ndarray     time vector of time domain signal of day of interest

    
    '''


    SR = int(trace[0].stats.sampling_rate)    # extracts samling rate to set decimate_factor

    ### as the amplitudes on a spectrogram calculated with mlab.specgram have a significant drop at the start/end, specgram calculations are made 
    ### for 26 hours. after the specgram is calculated, the central 24 hours of interest are cutted out to get daily spectrograms

    ## calculate values to trim for 26 hours of data
    time_win_to_use = trace[0].stats.starttime+24*3600                       # start middle day
    trim_start = time_win_to_use - (add_time)*3600                           # add one hour of day before
    trim_end = time_win_to_use + (24+add_time)*3600 - trace[0].stats.delta   # add one houre of day after

    ### preprocessing: 
    # - decimate from 500Hz to 100Hz
    # - filter between 5 to 49 Hz
    # - demean
    # - rotate data from DH1, DH2, DHZ to EHE, EHN, EHZ
    # - trim to get uniform time length for waveform data,
#    trace.decimate(factor=decimate_factor) 
#    trace.remove_response(inventory=inv, pre_filt=pre_filt, output="ACC",water_level=60, plot=False)
#    trace.filter('bandpass',freqmin=bandpass_fmin,freqmax=bandpass_fmax,zerophase=True)
#    trace.detrend(type='simple') # remove mean
#    trace.detrend(type='linear') # remove trend
    trace.trim(starttime=trim_start,endtime=trim_end,nearest_sample=True,pad=True, fill_value=0)

    ### hand over data from obpsy to numpy

    trace_data = np.array(trace[0].data)
    trace_data -= trace_data.mean()

    # hann (hanning) window is used for tapering the single time window for spectrum calculation in mlab.specgram
    specgram, specgram_freq_vec, specgram_time_vex = mlab.specgram(trace_data, Fs=SR, NFFT=nfft, noverlap=0)

    ff_value = np.zeros((np.shape(specgram)[1]))
    ff_amplitude =  np.zeros_like(ff_value)
    IH_max_amp = np.zeros_like(ff_value)
    trace_time_domain = np.zeros((len(trace_data)))

    ## calculation of fundamental frequency values by HarmonicProductSpectrum_algorithm, for futher information see ./func/def_master.py

    ff_value, ff_amplitude, IH_max_amp = hps_algo(input_data=specgram, freq_vector=specgram_freq_vec, nr_downsamp=nr_downsamp_HPS, min_freq=-1, max_freq=0.000001)

    ## saves time domain data of all component

    trace_time_domain= trace_data

    ff_amplitude = ff_amplitude[45:1125]
    ff_value = ff_value[45:1125]  
    IH_max_amp = IH_max_amp[45:1125]
    

    specgram = specgram[:,45:1125]
    specgram_time_vec = specgram_time_vex[45:1125]#-3600

    ## generate time line for plotting

    trace_data_one_day = trace_time_domain[3600*SR:len(trace_data)-3600*SR]
    trace_data_one_day_time_line = np.linspace(0,24,np.shape(trace_data_one_day)[0])

    ## Distinguish between signal/no signal for the plots of specgrams which are created later

    IHS = ff_amplitude/IH_max_amp

    ff_value_signal= ff_value.copy()
    ff_amplitude_signal= ff_amplitude.copy()
    ff_value_NOsignal= ff_value.copy()
    ff_amplitude_NOsignal= ff_amplitude.copy()

    for i in range(np.shape(ff_value)[0]):
        if IHS[i] < threshold_HSI:
            ff_value_signal[i] = np.nan
            ff_amplitude_signal[i] = np.nan
        if IHS[i] >= threshold_HSI:
            ff_value_NOsignal[i] = np.nan
            ff_amplitude_NOsignal[i] = np.nan
    return (ff_value_signal,ff_amplitude_signal,ff_value_NOsignal,ff_amplitude_NOsignal, specgram, specgram_time_vec ,specgram_freq_vec,trace_data_one_day,trace_data_one_day_time_line)

 
