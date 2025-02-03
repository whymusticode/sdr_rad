#!/usr/bin/env python

import h5py
import glob
import numpy
import matplotlib.pyplot as plt


fl = glob.glob("/data1/results/fm105.10/north/passivew/*/*.h5")
fl.sort()
#print(len(fl))

n_res = len(fl)-1
h=h5py.File(fl[0],"r")
lags_us = 1e6*h["lags"].value/100e3
range_km = numpy.copy(h["range"].value/1e3)
h.close()
M = numpy.zeros([len(range_km),n_res])
t = []
for i,f in enumerate(fl[0:n_res]):
#    print(".")
    h=h5py.File(f,"r")
    M[:,i]=numpy.copy(numpy.abs(h["lag_matrix"].value[:,0]))
    M[:,i]=M[:,i]/numpy.median(M[:,i])
    t.append(h["t"].value)
    h.close()
t=numpy.array(t)
#plt.pcolormesh(numpy.array(t),range_km,M)
plt.pcolormesh(t,range_km,10.0*numpy.log10(M),vmin=0.0)
plt.xlabel("Time (unix second)")
plt.ylabel("Range (km)")
plt.title("Passive radar zero lag estimate")
plt.xlim([numpy.min(t),numpy.max(t)])
plt.ylim([numpy.min(range_km),numpy.max(range_km)])
plt.colorbar()
plt.show()
