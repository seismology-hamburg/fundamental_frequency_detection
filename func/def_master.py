#!/bin/python3_6


import numpy as np
import numpy.ma as ma



def hps_algo(input_data, freq_range_input_data, min_freq, max_freq, nr_downsamp, amp_threshold=1e-45):
    '''HarmonicProductSpectrum Algorithum to obtain fundamental frequency  
    
    http://musicweb.ucsd.edu/~trsmyth/analysis/Harmonic_Product_Spectrum.html 
    
    INPUT
    
    input_data                          ndarray      spectrogram with harmonic signal on it
    freq_range_input_data               ndarray      frequency values calculated from mlab.specgram
    min_freq                            float        minimal frequency to mask frequency band, optional
    max_freq                            float        maximal frequency to mask frequency band
    nr_dwonsamp                         integer      value to downsamp 
    amp_samples_to_mute                 integer      ampltude values to mute at corner of every sepctrogram, if not set 8 is automatically choosen
    amp_threshold                       float        value to differ between harmonic signal and noise, optional   
    
    OUTPUT
    
    fund_freq_amplitude_all             ndarray     all amplitude values           
    fund_freq_amplitude                 ndarray     amplitude values over threshold
    fund_freq_value_all                 ndarray     all frequency values
    fund_freq_value                     ndarray     frequency values over threshold
    fund_freq_indi                      ndarray     all frequency value idices
    
    '''
    
    
    # mask array in given frequency band
    input_data_copy = input_data.copy()
    freq = freq_range_input_data
    
    lower_minval_freq = min_freq
    lower_maxval_freq = max_freq
    lower_freq_values = freq[np.where(np.logical_and(freq>lower_minval_freq, freq<lower_maxval_freq))]

    lower_first_indi_to_use = np.where(freq == lower_freq_values[0])[0]
    lower_last_indi_to_use = np.where(freq == lower_freq_values[-1])[0]

    input_data_masked = ma.array(input_data_copy)
    input_data_masked[int(lower_first_indi_to_use):int(lower_last_indi_to_use),:] = 0
    
    
    
    # initialize arrays to safe values
    # 2D arrays
    fund_freq_amp = np.zeros_like(input_data)
    harmonic_product_spectrum = np.zeros_like(input_data)
    #print(np.shape)
    
    # 1D arrays
    fund_freq_value = np.zeros((np.shape(input_data))[1])
    fund_freq_value_all = np.zeros_like(fund_freq_value)
    fund_freq_indi = np.zeros_like(fund_freq_value)
    fund_freq_amplitude = np.zeros_like(fund_freq_value)
    fund_freq_amplitude_all = np.zeros_like(fund_freq_value)
    fund_freq_amplitude_all_muted = np.zeros_like(fund_freq_value)

    
    spectrum_downsamp =  np.zeros(((np.shape(input_data_masked))[1],nr_downsamp, len(input_data_masked)))
    #print(np.shape(spectrum_downsamp), ' = spectrum_downsamp')
    
    
    # for loop over all columns
    for j in range((np.shape(input_data_masked))[1]):
        #print((np.shape(input_data_masked))[1])
        spectrum_signal = input_data_masked[:,j]
        #print(np.shape(spectrum_downsamp))
        spectrum_downsamp[j,0,:] = spectrum_signal
        #print(np.shape(spectrum_signal), ' = spectrum_signal')
    
        # for loop to downsample every column in spectrogram
        for i in range(nr_downsamp-1):
            nr_downsamp_array = np.linspace(2,nr_downsamp, nr_downsamp-1)
            #print(nr_downsamp_array)
            spectrum_downsamp[j,i+1,:][0:len(spectrum_signal[::int(nr_downsamp_array[i])])] = spectrum_signal[::int(nr_downsamp_array[i])]
            
            

        if nr_downsamp==4:
            harmonic_product_spectrum[:,j] = spectrum_downsamp[j,0,:] *  spectrum_downsamp[j,1,:] * spectrum_downsamp[j,2,:] * spectrum_downsamp[j,3,:]#* one_spec_downsamp[4]


        if nr_downsamp==3:
            harmonic_product_spectrum[:,j] = spectrum_downsamp[j,0,:] *  spectrum_downsamp[j,1,:] * spectrum_downsamp[j,2,:]# * spectrum_downsamp[j,3,:]#* one_spec_downsamp[4]

        if nr_downsamp==2:
            harmonic_product_spectrum[:,j] = spectrum_downsamp[j,0,:] *  spectrum_downsamp[j,1,:]# * spectrum_downsamp[j,2,:]# * spectrum_downsamp[j,3,:]#* one_spec_downsamp[4]
            
            
        fund_freq_value[j] = freq[np.argmax(harmonic_product_spectrum[:,j])]
        fund_freq_value_all[j] = freq[np.argmax(harmonic_product_spectrum[:,j])]
        #print(np.shape(harmonic_product_spectrum), ' = harmonic_product_spectrum')
        fund_freq_indi[j] = np.argmax(harmonic_product_spectrum[:,j])
        fund_freq_amplitude[j] = np.max(spectrum_signal[int(fund_freq_indi[j])])
        fund_freq_amplitude_all[j] = np.max(spectrum_signal[int(fund_freq_indi[j])])
        #fund_freq_amplitude_all_muted[j] = np.max(spectrum_signal[int(fund_freq_indi[j])])
        #print(10*np.log10(fund_freq_amplitude[::10]))
    
        # Amplitude Values for Threshold handling
    
        if fund_freq_amplitude[j] < amp_threshold:

            fund_freq_value[j] = np.nan
            #fund_freq_indi = np.nan
            fund_freq_amplitude[j] = np.nan    
            

            
    return fund_freq_amplitude_all, fund_freq_amplitude, fund_freq_value_all, fund_freq_value, fund_freq_indi 


def hps_algov3(input_data, freq_range_input_data, min_freq, max_freq, nr_downsamp, number_of_days, amp_samples_to_mute=12,amp_threshold=1e-45):
    '''HarmonicProductSpectrum Algorithum to obtain fundamental frequency  
    
    http://musicweb.ucsd.edu/~trsmyth/analysis/Harmonic_Product_Spectrum.html 
    
    INPUT
    
    input_data                          ndarray      spectrogram with harmonic signal on it
    freq_range_input_data               ndarray      frequency values calculated from mlab.specgram
    min_freq                            float        minimal frequency to mask frequency band, optional
    max_freq                            float        maximal frequency to mask frequency band
    nr_dwonsamp                         integer      value to downsamp 
    number_of_days                      integer      value of days appended together in spectrogram
    amp_samples_to_mute                 integer      ampltude values to mute at corner of every sepctrogram, if not set 8 is automatically choosen
    amp_threshold                       float        value to differ between harmonic signal and noise, optional   
    
    
    OUTPUT
    
    fund_freq_amplitude_all             ndarray     
    fund_freq_amplitude
    fund_freq_amplitude_all_muted
    fund_freq_indi
    fund_freq_value
    fund_freq_value_all
    
    '''
    
    
    # mask array in given frequency band
    input_data_copy = input_data.copy()
    freq = freq_range_input_data
    
    lower_minval_freq = min_freq
    lower_maxval_freq = max_freq
    lower_freq_values = freq[np.where(np.logical_and(freq>lower_minval_freq, freq<lower_maxval_freq))]

    lower_first_indi_to_use = np.where(freq == lower_freq_values[0])[0]
    lower_last_indi_to_use = np.where(freq == lower_freq_values[-1])[0]

    input_data_masked = ma.array(input_data_copy)
    input_data_masked[int(lower_first_indi_to_use):int(lower_last_indi_to_use),:] = 0
    
    
    
    # initialize arrays to safe values
    # 2D arrays
    fund_freq_amp = np.zeros_like(input_data)
    harmonic_product_spectrum = np.zeros_like(input_data)
    #print(np.shape)
    
    # 1D arrays
    fund_freq_value = np.zeros((np.shape(input_data))[1])
    fund_freq_value_all = np.zeros_like(fund_freq_value)
    fund_freq_indi = np.zeros_like(fund_freq_value)
    fund_freq_amplitude = np.zeros_like(fund_freq_value)
    fund_freq_amplitude_all = np.zeros_like(fund_freq_value)
    fund_freq_amplitude_all_muted = np.zeros_like(fund_freq_value)

    
    spectrum_downsamp =  np.zeros(((np.shape(input_data_masked))[1],nr_downsamp, len(input_data_masked)))
    #print(np.shape(spectrum_downsamp), ' = spectrum_downsamp')
    
    
    # for loop over all columns
    for j in range((np.shape(input_data_masked))[1]):
        #print((np.shape(input_data_masked))[1])
        spectrum_signal = input_data_masked[:,j]
        #print(np.shape(spectrum_downsamp))
        spectrum_downsamp[j,0,:] = spectrum_signal
        #print(np.shape(spectrum_signal), ' = spectrum_signal')
    
        # for loop to downsample every column in spectrogram
        for i in range(nr_downsamp-1):
            nr_downsamp_array = np.linspace(2,nr_downsamp, nr_downsamp-1)
            #print(nr_downsamp_array)
            spectrum_downsamp[j,i+1,:][0:len(spectrum_signal[::int(nr_downsamp_array[i])])] = spectrum_signal[::int(nr_downsamp_array[i])]
            
            

        if nr_downsamp==4:
            harmonic_product_spectrum[:,j] = spectrum_downsamp[j,0,:] *  spectrum_downsamp[j,1,:] * spectrum_downsamp[j,2,:] * spectrum_downsamp[j,3,:]#* one_spec_downsamp[4]


        if nr_downsamp==3:
            harmonic_product_spectrum[:,j] = spectrum_downsamp[j,0,:] *  spectrum_downsamp[j,1,:] * spectrum_downsamp[j,2,:]# * spectrum_downsamp[j,3,:]#* one_spec_downsamp[4]

        if nr_downsamp==2:
            harmonic_product_spectrum[:,j] = spectrum_downsamp[j,0,:] *  spectrum_downsamp[j,1,:]# * spectrum_downsamp[j,2,:]# * spectrum_downsamp[j,3,:]#* one_spec_downsamp[4]
            
            
        fund_freq_value[j] = freq[np.argmax(harmonic_product_spectrum[:,j])]
        fund_freq_value_all[j] = freq[np.argmax(harmonic_product_spectrum[:,j])]
        #print(np.shape(harmonic_product_spectrum), ' = harmonic_product_spectrum')
        fund_freq_indi[j] = np.argmax(harmonic_product_spectrum[:,j])
        fund_freq_amplitude[j] = np.max(spectrum_signal[int(fund_freq_indi[j])])
        fund_freq_amplitude_all[j] = np.max(spectrum_signal[int(fund_freq_indi[j])])
        fund_freq_amplitude_all_muted[j] = np.max(spectrum_signal[int(fund_freq_indi[j])])
        #print(10*np.log10(fund_freq_amplitude[::10]))
    
        # Amplitude Values for Threshold handling
    
        if fund_freq_amplitude[j] < amp_threshold:

            fund_freq_value[j] = np.nan
            #fund_freq_indi = np.nan
            fund_freq_amplitude[j] = np.nan
            
            
            
    for i in range(number_of_days):
        if i !=0:
            first_value = (int(len(fund_freq_amplitude_all_muted)/number_of_days*i)-amp_samples_to_mute)
            last_value = (int(len(fund_freq_amplitude_all_muted)/number_of_days*i)+amp_samples_to_mute+1)
            fund_freq_amplitude_all_muted[first_value:last_value] = np.nan
    fund_freq_amplitude_all_muted[:amp_samples_to_mute] = np.nan
    #fund_freq_amplitude_all_muted[:, int(np.shape(fund_freq_amplitude_all_muted)[1]/number_of_days)-amp_samples_to_mute:int(np.shape(fund_freq_amplitude_all_muted)[1]/number_of_days)+amp_samples_to_mute] = np.nan
    fund_freq_amplitude_all_muted[-amp_samples_to_mute:] = np.nan
    
            

            
    return fund_freq_amplitude_all, fund_freq_amplitude, fund_freq_amplitude_all_muted, fund_freq_indi, fund_freq_value, fund_freq_value_all 






def hps_algov2(input_data, freq_range_input_data, min_freq, max_freq, nr_downsamp, amp_threshold=1e-45):
    '''HarmonicProductSpectrum Algorithum to obtain fundamental frequency  
    
    http://musicweb.ucsd.edu/~trsmyth/analysis/Harmonic_Product_Spectrum.html 
    
    input_data              ndarray      spectrogram with harmonic signal on it
    freq_range_input_data   ndarray      frequency values calculated from mlab.specgram
    min_freq                float        minimal frequency to mask frequency band, optional
    max_freq                float        maximal frequency to mask frequency band
    nr_dwonsamp             integer      value to downsamp 
    amp_threshold           float        value to differ between harmonic signal and noise, optional   
    
    '''
    
    
    # mask array in given frequency band
    input_data_copy = input_data.copy()
    freq = freq_range_input_data
    
    lower_minval_freq = min_freq
    lower_maxval_freq = max_freq
    lower_freq_values = freq[np.where(np.logical_and(freq>lower_minval_freq, freq<lower_maxval_freq))]

    lower_first_indi_to_use = np.where(freq == lower_freq_values[0])[0]
    lower_last_indi_to_use = np.where(freq == lower_freq_values[-1])[0]

    input_data_masked = ma.array(input_data_copy)
    input_data_masked[int(lower_first_indi_to_use):int(lower_last_indi_to_use),:] = 0
    
    
    
    # initialize arrays to safe values
    # 2D arrays
    fund_freq_amp = np.zeros_like(input_data)
    harmonic_product_spectrum = np.zeros_like(input_data)
    #print(np.shape)
    
    # 1D arrays
    fund_freq_value = np.zeros((np.shape(input_data))[1])
    fund_freq_indi = np.zeros_like(fund_freq_value)
    fund_freq_amplitude = np.zeros_like(fund_freq_value)
    amp_fundamenteal_freq = np.zeros_like(fund_freq_value)

    
    spectrum_downsamp =  np.zeros(((np.shape(input_data_masked))[1],nr_downsamp, len(input_data_masked)))
    #print(np.shape(fund_freq_value), np.shape(fund_freq_indi), np.shape(amp_fundamenteal_freq))
    
    
    # for loop over all columns
    for j in range((np.shape(input_data_masked))[1]):
        #print((np.shape(input_data_masked))[1])
        spectrum_signal = input_data_masked[:,j]
        #print(np.shape(spectrum_downsamp))
        spectrum_downsamp[j,0,:] = spectrum_signal 
    
        # for loop to downsample every column in spectrogram
        for i in range(nr_downsamp-1):
            nr_downsamp_array = np.linspace(2,nr_downsamp, nr_downsamp-1)
            #print(nr_downsamp_array)
            spectrum_downsamp[j,i+1,:][0:len(spectrum_signal[::int(nr_downsamp_array[i])])] = spectrum_signal[::int(nr_downsamp_array[i])]
            
            

        if nr_downsamp==4:
            harmonic_product_spectrum[:,j] = spectrum_downsamp[j,0,:] *  spectrum_downsamp[j,1,:] * spectrum_downsamp[j,2,:] * spectrum_downsamp[j,3,:]#* one_spec_downsamp[4]


        if nr_downsamp==3:
            harmonic_product_spectrum[:,j] = spectrum_downsamp[j,0,:] *  spectrum_downsamp[j,1,:] * spectrum_downsamp[j,2,:]# * spectrum_downsamp[j,3,:]#* one_spec_downsamp[4]

        if nr_downsamp==2:
            harmonic_product_spectrum[:,j] = spectrum_downsamp[j,0,:] *  spectrum_downsamp[j,1,:]# * spectrum_downsamp[j,2,:]# * spectrum_downsamp[j,3,:]#* one_spec_downsamp[4]
            
            
        fund_freq_value[j] = freq[np.argmax(harmonic_product_spectrum[:,j])]
        fund_freq_indi[j] = np.argmax(harmonic_product_spectrum[:,j])
        fund_freq_amplitude[j] = np.max(harmonic_product_spectrum[:,j])
    
        # Amplitude Values for Threshold handling
    
        if fund_freq_amplitude[j] < amp_threshold:

            fund_freq_value[j] = np.nan
            #fund_freq_indi = np.nan
            fund_freq_amplitude[j] = np.nan

            
    return fund_freq_value, fund_freq_indi, fund_freq_amplitude




def hps_algov1(input_data, freq_range_input_data, min_freq, max_freq, nr_downsamp, amp_threshold=1e-45):
    '''HarmonicProductSpectrum Algorithum to obtain fundamental frequency  
    
    http://musicweb.ucsd.edu/~trsmyth/analysis/Harmonic_Product_Spectrum.html 
    
    input_data              ndarray      spectrogram with harmonic signal on it
    freq_range_input_data   ndarray      frequency values calculated from mlab.specgram
    min_freq                float        minimal frequency to mask frequency band, optional
    max_freq                float        maximal frequency to mask frequency band
    nr_dwonsamp             integer      value to downsamp 
    amp_threshold           float        value to differ between harmonic signal and noise, optional   
    
    '''
    
    
    # mask array in given frequency band
    input_data_copy = input_data.copy()
    freq = freq_range_input_data
    
    lower_minval_freq = min_freq
    lower_maxval_freq = max_freq
    lower_freq_values = freq[np.where(np.logical_and(freq>lower_minval_freq, freq<lower_maxval_freq))]

    lower_first_indi_to_use = np.where(freq == lower_freq_values[0])[0]
    lower_last_indi_to_use = np.where(freq == lower_freq_values[-1])[0]

    input_data_masked = ma.array(input_data_copy)
    input_data_masked[int(lower_first_indi_to_use):int(lower_last_indi_to_use),:] = 0
    
    
    
    # initialize arrays to safe values
    # 2D arrays
    fund_freq_amp = np.zeros_like(input_data)
    harmonic_product_spectrum = np.zeros_like(input_data)
    #print(np.shape)
    
    # 1D arrays
    fund_freq_value = np.zeros((np.shape(input_data))[1])
    fund_freq_indi = np.zeros_like(fund_freq_value)
    fund_freq_amplitude = np.zeros_like(fund_freq_value)
    amp_fundamenteal_freq = np.zeros_like(fund_freq_value)

    
    spectrum_downsamp =  np.zeros(((np.shape(input_data_masked))[1],nr_downsamp, len(input_data_masked)))
    #print(np.shape(fund_freq_value), np.shape(fund_freq_indi), np.shape(amp_fundamenteal_freq))
    
    
    # for loop over all columns
    for j in range((np.shape(input_data_masked))[1]):
        #print((np.shape(input_data_masked))[1])
        spectrum_signal = input_data_masked[:,j]
        #print(np.shape(spectrum_downsamp))
        spectrum_downsamp[j,0,:] = spectrum_signal 
    
        # for loop to downsample every column in spectrogram
        for i in range(nr_downsamp-1):
            nr_downsamp_array = np.linspace(2,nr_downsamp, nr_downsamp-1)
            #print(nr_downsamp_array)
            spectrum_downsamp[j,i+1,:][0:len(spectrum_signal[::int(nr_downsamp_array[i])])] = spectrum_signal[::int(nr_downsamp_array[i])]
        
        harmonic_product_spectrum[:,j] = spectrum_downsamp[j,0,:] *  spectrum_downsamp[j,1,:] * spectrum_downsamp[j,2,:] * spectrum_downsamp[j,3,:]#* one_spec_downsamp[4]

        fund_freq_value[j] = freq[np.argmax(harmonic_product_spectrum[:,j])]
        fund_freq_indi[j] = np.argmax(harmonic_product_spectrum[:,j])
        fund_freq_amplitude[j] = np.max(harmonic_product_spectrum[:,j])
    
        # Amplitude Values for Threshold handling
    
        if fund_freq_amplitude[j] < amp_threshold:

            fund_freq_value[j] = np.nan
            #fund_freq_indi = np.nan
            fund_freq_amplitude[j] = np.nan

            
    return fund_freq_value, fund_freq_indi, fund_freq_amplitude
