#!/usr/bin/env python
#
# Passive radar analysis
# incoherent integration using LPI
# (c) 2013, Juha Vierinen, MIT Haystack Observatory
#
import gdf
import scipy, scipy.constants, scipy.linalg, numpy
import h5py
import stuffr
import os, time, optparse, errno, shutil, glob
import matplotlib.pyplot as plt
from multiprocessing import Process
import pyfftw
import scipy.signal
import digital_rf_hdf5 as drf
import lp
import math
import sampler_util

def spectrogramc(x,window=1024,wf=scipy.signal.hamming):
    wfv = wf(window)
    Nwindow = int(math.floor(len(x)/window))
    res = numpy.zeros([Nwindow,window],dtype=numpy.complex64)
    for i in range(Nwindow):
        res[i,] = numpy.fft.fftshift(numpy.fft.fft(wfv*x[i*window + numpy.arange(window)]))
    return(res)

def find_bins(cfreqs,freqs):
    bins = []
    for f in freqs:
        i = numpy.argmin(numpy.abs(cfreqs - f))
        bins.append(i)
        print(numpy.abs(cfreqs[i] - f))
    return(numpy.array(bins,dtype=numpy.int))

# polyphase filterbank window
def get_window(n=100,p=4):
    w=numpy.zeros([n*p],dtype=numpy.complex64)
    w[0:(p/2)]=1.0
    w[(n*p-(p/2)+1):(n*p)]=1.0
    h = numpy.array(scipy.signal.chebwin(n*p,100.0)*numpy.fft.fftshift(numpy.fft.ifft(w)),dtype=numpy.complex64)
    #h = numpy.array(numpy.fft.fftshift(numpy.fft.ifft(w)),dtype=numpy.complex64)
    return(h)
#    numpy.fft.ifft(w)


d = drf.read_hdf5("/data1/hay_n200_passive_radar")
sr = long(d.get_metadata("fm_south")["sample_rate"].value)
cf = d.get_metadata("fm_south")["center_frequencies"].value[0]
t0 = long(d.get_metadata("fm_south")["t0"].value)*100000

#dcfreq = 101e6-105.7e6
freqs = numpy.array([105.1,104.9,104.5,104.1,103.3,102.5,100.3,100.1,101.1,97.7,97.3])*1e6
#freqs = numpy.array([101.1,100.9,105.7,105.1])*1e6#,104.9,104.5,104.1,103.3,102.5,100.3,100.1,97.7,97.3])*1e6

dec = 100  # 100 kHz out
pl = 8     # 4 window pfb
window = get_window(n=dec,p=pl)

dcfreq=0.0
sr2 = sr/dec
cfreqs = numpy.fft.fftshift((numpy.arange(dec)-dec/2)*sr2 + cf)
bins = find_bins(cfreqs,freqs)
labels = ["%1.2f"%(f/1e6) for f in freqs]
print(cfreqs)
print(labels)
print(bins)

plot_passband=True
if plot_passband:
    plt.plot(numpy.linspace(-sr/2.0,sr/2.0,num=dec*pl)/1e3,10.0*numpy.log10(numpy.abs(numpy.fft.fftshift(numpy.fft.fft(window)))**2.0))
    plt.xlabel("kHz")
    plt.ylabel("dB")
    plt.show()



# downconversion vector
ws = 10000000
zexp=numpy.array(numpy.exp(1.0j*2.0*numpy.pi*dcfreq*numpy.arange(ws+pl*dec,dtype=numpy.float64)/10e6),dtype=numpy.complex64)
zstep=numpy.array(numpy.exp(1.0j*2.0*numpy.pi*dcfreq*(ws+pl*dec + 1.0)/10e6),dtype=numpy.complex64)
phase=zexp[0]
b = d.get_bounds("fm_south")
nwin = long(numpy.floor((b[1]-b[0])/ws)-1)

zn_lp = numpy.zeros([len(bins),ws/dec],dtype=numpy.complex64)
zs_lp = numpy.zeros([len(bins),ws/dec],dtype=numpy.complex64)

idx = b[0]
debug_spec=False
debug_plot=False

if not debug_plot:
    dir_s=[]
    dir_n=[]
    wr_s=[]
    wr_n=[]
    for i,f in enumerate(freqs):
        dir_s.append("/data1/results/fm%s/south"%(labels[i]))
        dir_n.append("/data1/results/fm%s/north"%(labels[i]))
        os.system("rm -Rf %s"%(dir_s[i]))
        os.system("rm -Rf %s"%(dir_n[i]))
        os.system("mkdir -p %s"%(dir_s[i]))
        os.system("mkdir -p %s"%(dir_n[i]))
        wr_s.append(drf.write_hdf5_channel(dir_s[i],"f",int(sr/dec),3600,t0,int(sr/dec),"ASDSAD"))
        wr_n.append(drf.write_hdf5_channel(dir_n[i],"f",int(sr/dec),3600,t0,int(sr/dec),"ASDSAD"))
        sampler_util.write_metadata_drf(dir_s[i],1,[f],t0/100e3,dtype="<f4",itemsize=8,sr=sr/dec,extra_keys=["pfb"],extra_values=["1.0"])
        sampler_util.write_metadata_drf(dir_n[i],1,[f],t0/100e3,dtype="<f4",itemsize=8,sr=sr/dec,extra_keys=["pfb"],extra_values=["1.0"])

for i in range(nwin):
    if debug_spec:
        z=d.read_vector_c81d(idx,ws,"fm_south")
        S=stuffr.spectrogram(z,window=10000)
        plt.subplot(211)
        plt.pcolormesh(numpy.linspace(101e6-5e6,101e6+5e6,num=10000)/1e6,numpy.arange(1000),10.0*numpy.log10(S))
        plt.colorbar()

        z=d.read_vector_c81d(idx,ws,"fm_north")
        S=stuffr.spectrogram(z,window=10000)
        plt.subplot(212)
        plt.pcolormesh(numpy.linspace(101e6-5e6,101e6+5e6,num=10000)/1e6,numpy.arange(1000),10.0*numpy.log10(S))
        plt.colorbar()
        plt.show()

    z=phase*zexp*d.read_vector_c81d(idx,ws+pl*dec,"fm_south")

    a = lp.pfb(z, ws, zs_lp, dec, window, pl, bins)

    z=phase*zexp*d.read_vector_c81d(idx,ws+pl*dec,"fm_north")

    a = lp.pfb(z, ws, zn_lp, dec, window, pl, bins)

    if not debug_plot:
        for i,b in enumerate(bins):
            data_f = numpy.zeros((ws/dec, 2), dtype="f")
            data_f[:,0]=zn_lp.real[i,0:(ws/dec)]
            data_f[:,1]=zn_lp.imag[i,0:(ws/dec)]
            wr_n[i].rf_write(data_f)
            
            data_f[:,0]=zs_lp.real[i,0:(ws/dec)]
            data_f[:,1]=zs_lp.imag[i,0:(ws/dec)]
            wr_s[i].rf_write(data_f)


#class write_hdf5_channel:
#    """The class write_hdf5_channel is an object used to write rf data to Hdf5 files as specified
#    in the http://www.haystack.mit.edu/pipermail/rapid-dev/2014-February/000273.html email thread.
#    """
#    
#    def __init__(self, directory, dtype_str, samples_per_file, files_per_directory, start_global_index, sample_rate, uuid_str,
#                 compression_level=0, checksum=False, is_complex=True, num_subchannels=1, marching_periods=True):

    if debug_plot:
        for i,b in enumerate(bins):
            SS=spectrogramc(zs_lp[i,:])
            SN=spectrogramc(zn_lp[i,:])
        
            plt.subplot(121)
            plt.pcolormesh(10.0*numpy.log10(numpy.abs(SS)**2.0))
            #plt.pcolormesh(numpy.abs(SS*numpy.conj(SN)))
            plt.title("%s"%(labels[i]))
            plt.colorbar()
        
            plt.subplot(122)
            plt.pcolormesh(10.0*numpy.log10(numpy.abs(SN)**2.0))
            #plt.pcolormesh(numpy.angle(SS*numpy.conj(SN)))
            plt.colorbar()
            plt.show()
    phase=phase*zstep
    idx+=ws
    
    

#    
#    ws = drf.write_hdf5_channel("/data1/results/fm_south%s"%(label),"f",100000,3600,t0)

#class write_hdf5_channel:
#    """The class write_hdf5_channel is an object used to write rf data to Hdf5 files as specified
#    in the http://www.haystack.mit.edu/pipermail/rapid-dev/2014-February/000273.html email thread.
#    """
#    
#    def __init__(self, directory, dtype_str, samples_per_file, files_per_directory, start_global_index, sample_rate, uuid_str,
#                 compression_level=0, checksum=False, is_complex=True, num_subchannels=1, marching_periods=True):
