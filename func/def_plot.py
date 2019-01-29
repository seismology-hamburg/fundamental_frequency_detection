#!/bin/python3_6


import numpy as np
import matplotlib.pyplot as plt


def plot_data(components, trace, trace_time_vector,  specgram_time_vector, specgram_freq_vector,  specgram, freq_values, amplitude_values, path, station, year, day):
    #i = 0
    #while i < np.shape(specgram)[0]:
    gridsize = (55, 9)
    fig = plt.figure(figsize=(32, 22))
    fig.subplots_adjust(hspace=0, wspace=0)
        #fig.suptitle(station + '.' + 'starttime: ' +  first_date)
        #fig.suptitle(station + '.' + component + '_' + 'starttime: ' +  str(trace[0].stats.starttime)[0:19])

    ## v_component
    ax1 = plt.subplot2grid(gridsize, (0, 0), colspan=7, rowspan=2)
    ax2 = plt.subplot2grid(gridsize, (2, 0), colspan=7, rowspan=8, sharex=ax1)
    ax3 = plt.subplot2grid(gridsize, (10, 0), colspan=7, rowspan=3)
        
        
    ## h_component1

    ax4 = plt.subplot2grid(gridsize, (14, 0), colspan=7, rowspan=2, sharey=ax1)
    ax5 = plt.subplot2grid(gridsize, (16, 0), colspan=7, rowspan=8, sharex=ax1)
    ax6 = plt.subplot2grid(gridsize, (24, 0), colspan=7, rowspan=3, sharex=ax3, sharey=ax3)
        
    ## h_component2

    ax7 = plt.subplot2grid(gridsize, (28, 0), colspan=7, rowspan=2, sharey=ax1)
    ax8 = plt.subplot2grid(gridsize, (30, 0), colspan=7, rowspan=8, sharex=ax1)
    ax9 = plt.subplot2grid(gridsize, (38, 0), colspan=7, rowspan=3, sharex=ax3, sharey=ax3)
        

    ## hydrophone

    ax10 = plt.subplot2grid(gridsize, (42, 0), colspan=7, rowspan=2)
    ax11 = plt.subplot2grid(gridsize, (44, 0), colspan=7, rowspan=8, sharex=ax1)
    ax12 = plt.subplot2grid(gridsize, (52, 0), colspan=7, rowspan=3, sharex=ax3, sharey=ax3)


    ax1.plot(trace_time_vector,trace[0,:], color='maroon')
    ax1.set_xlim(0,24)
    ax1.title.set_text(components[0])
    ax1.set_ylim(-.005,.005)
    ax1.set_ylabel('Acceleration')
    ax1.set_xticks([])


    ax2.pcolormesh(specgram_time_vector/3600, specgram_freq_vector, 10*np.log10(specgram[0,:,:]))
    ax2.plot(specgram_time_vector/3600,freq_values[0,:],',', color='dodgerblue')
    #ax1.plot(time_stack/3600/24,freq_values_all_appended[0,:],',', color='dodgerblue')
    #ax2.set_ylim(0,5.5)
    ax2.set_ylabel('Frequency [Hz]')
    #ax2.set_xlim(0,24)
    #ax2.title.set_text(components[0])
    ax2.set_xticks([]) 


    #ax3.plot(specgram_time_vector/3600, 10*np.log10(amplitude_values_all[0,:]),',', color='red')
    #ax3.plot(time_stack/3600/24,10*np.log10(amplitude_values_all_appended[0,:]),',', color='red')
    ax3.plot(specgram_time_vector/3600,10*np.log10(amplitude_values[0,:]),',', color='maroon')
    ax3.axhline(y=-100, color='red', linewidth=.6, label='-100 dB')
    mean0 = str(10*np.log10(np.nanmean(amplitude_values[0,:])))
    ax3.axhline(y=10*np.log10(np.nanmean(amplitude_values[0,:])), color='lime', lw=.3, label=mean0[0:4])
    ax3.set_ylabel('Amplitude [dB]')
    ax3.set_ylim(-130,-60)
    ax3.set_xlim(0, 24)
    ax3.legend(bbox_to_anchor=(1.01, .6), loc=2, borderaxespad=0.)


    ax4.plot(trace_time_vector,trace[1,:], color='maroon')
    ax4.set_xlim(0,24)
    ax4.title.set_text(components[1])
    ax4.set_ylabel('Acceleration')
    ax4.set_xticks([])


    ax5.pcolormesh(specgram_time_vector/3600, specgram_freq_vector, 10*np.log10(specgram[1,:,:]))
    ax5.plot(specgram_time_vector/3600,freq_values[1,:],',', color='dodgerblue')
    #ax1.plot(time_stack/3600/24,freq_values_all_appended[0,:],',', color='dodgerblue')
    #ax2.set_ylim(0,5.5)
    ax5.set_ylabel('Frequency [Hz]')
    #a2.set_xlim(0,24)
    #ax2.title.set_text(components[0])
    ax5.set_xticks([]) 


    #ax3.plot(specgram_time_vector/3600, 10*np.log10(amplitude_values_all[0,:]),',', color='red')
    #ax3.plot(time_stack/3600/24,10*np.log10(amplitude_values_all_appended[0,:]),',', color='red')
    ax6.plot(specgram_time_vector/3600,10*np.log10(amplitude_values[0,:]),',', color='maroon')
    #a3.axhline(y=10*np.log10(threshold), color='red', linewidth=.6, label=10*np.log10(threshold))
    ax6.axhline(y=-100, color='red', linewidth=.6, label='-100 dB')
    mean0 = str(10*np.log10(np.nanmean(amplitude_values[1,:])))
    ax6.axhline(y=10*np.log10(np.nanmean(amplitude_values[1,:])), color='lime', lw=.3, label=mean0[0:4])
    ax6.set_ylabel('Amplitude [dB]')
    ax6.set_xlim(0, 24)
    ax6.legend(bbox_to_anchor=(1.01, .6), loc=2, borderaxespad=0.)

    ax7.plot(trace_time_vector,trace[2,:], color='maroon')
    ax7.set_xlim(0,24)
    ax7.title.set_text(components[2])
    ax7.set_ylabel('Acceleration')
    ax7.set_xticks([])


    ax8.pcolormesh(specgram_time_vector/3600, specgram_freq_vector, 10*np.log10(specgram[2,:,:]))
    ax8.plot(specgram_time_vector/3600,freq_values[2,:],',', color='dodgerblue')
    #ax1.plot(time_stack/3600/24,freq_values_all_appended[0,:],',', color='dodgerblue')
    #ax2.set_ylim(0,5.5)
    ax8.set_ylabel('Frequency [Hz]')
    #a2.set_xlim(0,24)
    #ax2.title.set_text(components[0])
    ax8.set_xticks([]) 

    #ax3.plot(specgram_time_vector/3600, 10*np.log10(amplitude_values_all[0,:]),',', color='red')
    #ax3.plot(time_stack/3600/24,10*np.log10(amplitude_values_all_appended[0,:]),',', color='red')
    ax9.plot(specgram_time_vector/3600,10*np.log10(amplitude_values[2,:]),',', color='maroon')
    ax9.axhline(y=-100, color='red', linewidth=.6, label='-100 dB')
    #a3.axhline(y=10*np.log10(threshold), color='red', linewidth=.6, label=10*np.log10(threshold))
    mean0 = str(10*np.log10(np.nanmean(amplitude_values[2,:])))
    ax9.axhline(y=10*np.log10(np.nanmean(amplitude_values[2,:])), color='lime', lw=.3, label=mean0[0:4])
    ax9.set_ylabel('Amplitude [dB]')
    ax9.set_xlim(0, 24)
    ax9.legend(bbox_to_anchor=(1.01, .6), loc=2, borderaxespad=0.)

    ax10.plot(trace_time_vector,trace[3,:], color='maroon')
    ax10.set_xlim(0,24)
    ax10.set_ylim(0.5,0.5)
    ax10.title.set_text(components[3])
    ax10.set_ylabel('Acceleration')
    ax10.set_xticks([])


    ax11.pcolormesh(specgram_time_vector/3600, specgram_freq_vector, 10*np.log10(specgram[3,:,:]))
    ax11.plot(specgram_time_vector/3600,freq_values[3,:],',', color='dodgerblue')
    #ax1.plot(time_stack/3600/24,freq_values_all_appended[0,:],',', color='dodgerblue')
    #ax2.set_ylim(0,5.5)
    ax11.set_ylabel('Frequency [Hz]')
    #a2.set_xlim(0,24)
    #ax2.title.set_text(components[0])
    ax11.set_xticks([]) 

    #ax3.plot(specgram_time_vector/3600, 10*np.log10(amplitude_values_all[0,:]),',', color='red')
    #ax3.plot(time_stack/3600/24,10*np.log10(amplitude_values_all_appended[0,:]),',', color='red')
    ax12.plot(specgram_time_vector/3600,10*np.log10(amplitude_values[3,:]),',', color='maroon')
    ax12.axhline(y=-100, color='red', linewidth=.6, label='-100 dB')
    #a3.axhline(y=10*np.log10(threshold), color='red', linewidth=.6, label=10*np.log10(threshold))
    mean0 = str(10*np.log10(np.nanmean(amplitude_values[3,:])))
    ax12.axhline(y=10*np.log10(np.nanmean(amplitude_values[3,:])), color='lime', lw=.3, label=mean0[0:4])
    ax12.set_ylabel('Amplitude [dB]')
    ax12.set_xlim(0, 24)
    ax12.legend(bbox_to_anchor=(1.01, .6), loc=2, borderaxespad=0.)

        #i += 1
    plt.suptitle((station+'.'+str(year)+' Day Nr: '+str(day)), x=.45, y=.91)    
    plt.savefig(path+'/'+station+'/1B.'+station+'.'+str(year)+'.'+day+'.png')
    plt.close('all')
    fig.clf()
