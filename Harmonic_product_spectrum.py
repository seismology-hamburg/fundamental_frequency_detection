
# coding: utf-8

import numpy as np
from obspy import read, read_inventory
#from obspy.io.xseed import Parser
from obspy.imaging.spectrogram import _nearest_pow_2
from scipy import signal
import obspy
import matplotlib.pyplot as plt
from IPython.display import Image

from matplotlib import mlab, transforms
plt.rcParams['image.cmap'] = 'inferno'
plt.rcParams['agg.path.chunksize'] = 10000
import math as M
import matplotlib as mpl
import numpy.ma as ma
import time as time
import os, gc, re

from func.def_master import hps_algo
from func.def_plot import plot_data



station_list = ['KNR16', 'KNR22']
#components = ['BHZ']

    
#components = ['BHZ', 'BH1', 'BH2', 'BDH']## Component
#components = ['BHZ']
year = 2016
#day = '065'
startday = 290

#pre_filt=(.001,.005,12,15) ## Parameter for filtering during response remove
pre_filt=(.05,.009,10,12) ## Parameter for filtering during response remove
#data_decimate_factor = 2
freq_min_bandpass = .5 ## filtercorners for filtering after response remove
freq_max_bandpass = 10 ## filtercorners for filtering after response remove
#per_lap = .9
#path_npy_files = '/home/david/AWI/Results/'+str(year)+'/'+str(freq_min_bandpass) +'Hz/'
#path_plots = path_npy_files

## Test environment, needs to be deleted...

path_npy_files = '/home/david/AWI/Results/TEST/'+str(year)+'/'
path_plots = path_npy_files

if os.path.isdir(path_npy_files) == False:
    os.mkdir(path_npy_files)
    
    
    
    
## original to calculate with gaps inbetween days
#start_day = str(day)


number_of_days = 2

threshold = -11 ## to delete amp values under certain amplitude threshold, 
threshold = 10**(threshold) ## in hps_algo; as I save all values its obsolete

save_amp_f=True #saves amplitude vaules, frequency values as .npy files in to path_npy_files
save_time_spec=False #saves spectrogram, freq and time vector of specgram and time_domain  as .npy files in to path_npy_files (disk space consuming!!!)
plot_it=True # to plot specgram with picked values for freq and amp in .png file

add_time = 1 # hours to add before/after day in hour
#per_lap = .1




#freq_minimal_value = 3
#distance = 120
#maxima= 1
#decimate_factor = 1e-15



for l in range(len(station_list)):
    day = startday
    station = station_list[l]
    station_value = re.findall('\d*\.?\d+',station)
    if int(station_value[0]) <=  17 or int(station_value[0]) == 19:
        components = ['BHZ', 'BH1', 'BH2', 'BDH']## Component
    
    else:
        components = ['HHZ', 'HH1', 'HH2', 'HDH']## Component
        
        
    if os.path.isdir(path_npy_files+station) == False:
        os.mkdir(path_npy_files+station)
        
    for i in range(number_of_days):
        for j in range(len(components)):
            print(day)
            ## Depends on workstation 
            #trace = read('/home/david/AWI/'+year+'/1B/' + station + '/' + components[j] + '.D/1B.' + station + '..' + components[j] + '.D.' + year + '.' + str(int(day)))
            #trace += read('/home/david/AWI/'+year+'/1B/' + station + '/' + components[j] + '.D/1B.' + station + '..' + components[j] + '.D.' + year + '.' + str(int(day)+1))
            #trace += read('/home/david/AWI/'+year+'/1B/' + station + '/' + components[j] + '.D/1B.' + station + '..' + components[j] + '.D.' + year + '.' + str(int(day)-1))
            trace = read('/home/david/AWI/'+str(year)+'/1B/' + station + '/' + components[j] + '.D/1B.' + station + '..' + components[j] + '.D.' + str(year) + '.' + '{:03}'.format(int(day)))
            if day==366:
                trace += read('/home/david/AWI/'+str(year)+'/1B/' + station + '/' + components[j] + '.D/1B.' + station + '..' + components[j] + '.D.' + str(year) + '.' + '{:03}'.format(int(day)-1))
                trace += read('/home/david/AWI/'+str(year+1)+'/1B/' + station + '/' + components[j] + '.D/1B.' + station + '..' + components[j] + '.D.' + str(year+1) + '.' + '001')
            elif day==1:
                trace += read('/home/david/AWI/'+str(year-1)+'/1B/' + station + '/' + components[j] + '.D/1B.' + station + '..' + components[j] + '.D.' + str(year-1) + '.' + '366')
                trace += read('/home/david/AWI/'+str(year)+'/1B/' + station + '/' + components[j] + '.D/1B.' + station + '..' + components[j] + '.D.' + str(year) + '.' + '{:03}'.format(int(day)+1))
            else:
                trace += read('/home/david/AWI/'+str(year)+'/1B/' + station + '/' + components[j] + '.D/1B.' + station + '..' + components[j] + '.D.' + str(year) + '.' + '{:03}'.format(int(day)+1))
                trace += read('/home/david/AWI/'+str(year)+'/1B/' + station + '/' + components[j] + '.D/1B.' + station + '..' + components[j] + '.D.' + str(year) + '.' + '{:03}'.format(int(day)-1))


            samp_rate = trace[0].stats.sampling_rate
            print(samp_rate)

            if samp_rate==50:
                data_decimate_factor = 2
                print('decimation with factor of 2')
            else:
                data_decimate_factor = 4
                print('decimation with factor of 4')
            

            # read in data for day before/after day of interest
            #trace = read('/data/cen/u254/Essing/AWI/KNIPA_NEW/' + year + '/1B/' + station + '/' + components[j] + '.D/1B.' + station + '..' + components[j] + '.D.' + year + '.' + str(int(day)))
            #trace += read('/data/cen/u254/Essing/AWI/KNIPA_NEW/' + year + '/1B/' + station + '/' + components[j] + '.D/1B.' + station + '..' + components[j] + '.D.' + year + '.' + str(int(day)+1))
            #trace += read('/data/cen/u254/Essing/AWI/KNIPA_NEW/' + year + '/1B/' + station + '/' + components[j] + '.D/1B.' + station + '..' + components[j] + '.D.' + year + '.' + str(int(day)-1))
            trace.merge(method=1)
            #print(trace[0].stats)

            ## calculate values to trim
            time_win_to_use = trace[0].stats.starttime+24*3600
            trim_start = time_win_to_use - (add_time)*3600
            trim_end = time_win_to_use + (24+add_time)*3600 - trace[0].stats.delta
            #print(trim_start, trim_end)

            trace.trim(trim_start, trim_end)
            #print(trace[0].stats)

            inv = read_inventory('/home/david/AWI/KNIPAS_dataless/RESP.1B.' + station + '..' + components[j])
            #resp_file = Parser('/home/david/AWI/KNIPAS_dataless/RESP.1B.' + station + '..' + components[j])
            #resp_file = ('/data/cen/u254/Essing/AWI/Dataless_SEEDS/RESP.1B.' + station + '..' + components[j])
            trace = trace.decimate(factor=data_decimate_factor)
            fs = int(trace[0].stats.sampling_rate)
            trace_simu = trace.copy()
            trace_simu.remove_response(inventory=inv, pre_filt=pre_filt, output="ACC",water_level=60, plot=False)
            #trace_simu.simulate(paz_remove=None, pre_filt=pre_filt, seedresp={'filename': resp_file, 'units': 'ACC'})
            trace_simu_filt = trace_simu.copy()
            trace_simu_filt = trace_simu_filt.filter('bandpass',freqmin=freq_min_bandpass, freqmax=freq_max_bandpass)
            trace_data = np.array(trace_simu_filt[0].data)
            trace_data -= trace_data.mean()




            #npts = len(trace_data)
            #wlen = 80 # in seconds
            #nfft = int(_nearest_pow_2(wlen * fs)) # number of datapoints used in each block for fft
            #nlap = int(nfft * float(per_lap))
            nfft = 2000
            pad_to = 2048

            # hanning window is set to default in mlab.specgram
            specgram, specgram_freq, specgram_time = mlab.specgram(trace_data, Fs=fs, NFFT=nfft, noverlap=0)
            #print(np.shape(specgram))
            #specgram_BHZ_masekd = ma.array(trace_data)
            #specgram_BHZ_masekd[0:1000,:] = ma.masked
            #specgram_BHZ_masekd[1300:np.shape(specgram_BHZ)[0],:] = ma.masked

            ### time arrangements

            if j==0:
                #specgram_appended_all_comps =  np.zeros((len(components), np.shape(specgram)[0], np.shape(specgram)[1]*number_of_days))
                first_date_1 = trace[0].stats.endtime.ctime()
                first_date = str(first_date_1[-4:]), '-', (first_date_1[4:7]), str(first_date_1[8:10])
                #all_days_amp_values = np.zeros((number_of_days, len(amplitude_values)))
                #all_days_freq_values = np.zeros((number_of_days, len(frequency_values)))
                time_stack = specgram_time


                fund_freq_value = np.zeros((len(components), np.shape(specgram)[1]))
                fund_freq_value_all =  np.zeros_like(fund_freq_value)
                fund_freq_indi = np.zeros_like(fund_freq_value)
                fund_freq_amplitude =  np.zeros_like(fund_freq_value)
                fund_freq_amplitude_all = np.zeros_like(fund_freq_value)
                trace_time_domain = np.zeros((len(components), len(trace_data)))
                #fund_freq_amplitude_all_muted = np.zeros_like(fund_freq_value_appended)


            #if i!=0 and j==0:
            #    time_stack = np.append(time_stack, time_stack[-1]+specgram_time)        



            #if i==0:      
            #    specgram_appended = specgram        
            #if i!=0:
            #    specgram_appended = np.append(specgram_appended, specgram, axis=1)

            if j==0:
                specgram_all_comps = np.zeros((len(components), np.shape(specgram)[0], np.shape(specgram)[1]))
                #specgram_appended = 0
                #print('I was here')
            specgram_all_comps[j,:,:] = specgram


            #print(np.shape(specgram_appended))

            #print(np.shape(fund_freq_value_appended), np.shape(fund_freq_indi_appended), np.shape(fund_freq_amplitude_appended))





            fund_freq_amplitude_all[j,:], fund_freq_amplitude[j,:], fund_freq_value_all[j,:], fund_freq_value[j,:], fund_freq_indi[j,:] = hps_algo(input_data=specgram, freq_range_input_data=specgram_freq, min_freq=-1, max_freq=0.000001, nr_downsamp=4, amp_threshold=threshold)
            trace_time_domain[j,:] = trace_data

        fund_freq_amplitude_all = fund_freq_amplitude_all[:,45:1125]
        fund_freq_amplitude = fund_freq_amplitude[:,45:1125]
        #fund_freq_amplitude_all_muted = fund_freq_amplitude_all_muted[:,45:1125]
        fund_freq_indi = fund_freq_indi[:,45:1125]
        fund_freq_value = fund_freq_value[:,45:1125]
        fund_freq_value_all = fund_freq_value_all[:,45:1125]

        specgram_all_comps = specgram_all_comps[:,:,45:1125]

        trace_data_one_day = trace_time_domain[:,3600*fs:len(trace_data)-3600*fs]
        trace_data_one_day_time_line = np.linspace(0,24,np.shape(trace_data_one_day)[1])


        if save_amp_f==True:
            np.save(file=(path_npy_files+'/'+station+'/1B.'+station+'.'+str(year)+'.'+str(day)+'_amp_values'), arr=fund_freq_amplitude_all)
            np.save(file=(path_npy_files+'/'+station+'/1B.'+station+'.'+str(year)+'.'+str(day)+'_freq_values'), arr=fund_freq_value_all)
            
        if save_time_spec==True:
            np.save(file=(path_npy_files+'/'+station+'/1B.'+station+'.'+str(year)+'.'+str(day)+'_specgrams'), arr=specgram_all_comps)
            np.save(file=(path_npy_files+'/'+station+'/1B.'+station+'.'+str(year)+'.'+str(day)+'_trace'), arr=trace_data_one_day)
            #fund_freq_amplitude_all_appended[j,:],  fund_freq_amplitude_appended[j,:],fund_freq_amplitude_all_muted_appended[j,:], fund_freq_indi_appended[j,:],  fund_freq_value_appended[j,:],  fund_freq_value_all_appended[j,:] = hps_algo(input_data=specgram_appended, freq_range_input_data=specgram_freq, min_freq=-1, max_freq=0.0001, nr_downsamp=4,number_of_days=number_of_days, amp_threshold=threshold)
            #fund_freq_value_appended[j,:], fund_freq_indi_append ed[j,:], fund_freq_amplitude_appended[j,:]  = hps_algo(input_data=specgram_appended, freq_range_input_data=freq, min_freq=-1, max_freq=0.0001, nr_downsamp=4, amp_threshold=threshold)

        middle_day_specgram_time = specgram_time[45:1125]-3600

        if plot_it==True:
            plot_data(components=components, trace=trace_data_one_day, trace_time_vector=trace_data_one_day_time_line,  specgram_time_vector=middle_day_specgram_time, specgram_freq_vector=specgram_freq,  specgram=specgram_all_comps, freq_values=fund_freq_value_all, amplitude_values=fund_freq_amplitude_all, path=path_plots, station=station, year=year,day=str(day))
        plt.close('all')
        gc.collect()
        day = day + 1
        #day = (("%03d" % day)) 
        ### Get amplitude at points in specgram without downsampling



    if save_time_spec==True:
        np.save(file=(path_npy_files+'/'+station+'/1B.'+station+'._time'), arr=middle_day_specgram_time)    
        np.save(file=(path_npy_files+'/'+station+'/1B.'+station+'._freq'), arr=specgram_freq)    
        np.save(file=(path_npy_files+'/'+station+'/1B.'+station+'._trace_time_line'), arr=trace_data_one_day_time_line)

