#!/usr/bin/python
#
# plotting routines for passive radar results
#
import matplotlib
matplotlib.use("Agg")

import itertools
import stuffr
import glob, os, numpy, h5py
import matplotlib.pyplot as plt
import matplotlib.colors
import os, time, optparse, errno, math

def get_combinations(n_channels):
    cs_iter = itertools.combinations(numpy.arange(n_channels),2)
    n_baselines=n_channels*(n_channels-1)/2
    cc_idx = numpy.zeros([n_baselines,2],dtype=numpy.int)
    for i,cs_pair in enumerate(cs_iter):
        cc_idx[i,:]=cs_pair
    return(cc_idx)

# mkdir -p sucks with python (this was borrowed from stackoverflow.com)
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def phase2imgcolors(s):
    L = s.shape[0]
    N = s.shape[1]
    hsv = numpy.zeros([L,N,3])

    # hue
    hsv[:,:,0] = (numpy.angle(s)+numpy.pi)/(2.0*numpy.pi)
    # sat (todo: calc spectral width)
    ass = stuffr.comprz_dB(numpy.abs(s)**2.0,fr=0.1)
    hsv[:,:,1] = 1.0
    # lum
    hsv[:,:,2] = (ass-numpy.min(ass))/(numpy.max(ass)-numpy.min(ass))
#    idx = numpy.where(hsv[:,:,2] < 0.4)
#    hsv[idx[0],idx[1],1] = 0.1

    rgb_image = matplotlib.colors.hsv_to_rgb(hsv)
    return(rgb_image)

def normalize(x):
    return((x-numpy.min(x))/(numpy.max(x)-numpy.min(x)))

# plot all individual spectra in directory
def plot_all_passive(dirn,prefix="passive",pprefix="passive/plot",n_incoh_int=10,dB_max=30.0, quiet=False, cc0=0,cc1=1):
    fl = glob.glob("%s/%s/cohin*.hdf5"%(dirn,prefix))
    r0 = h5py.File(fl[0],'r')
    CS = numpy.copy( r0["spec"][0,:,:] )
    CS[:,:]=0.0
    SH = []
    os.system("mkdir -p %s/%s"%(dirn,pprefix))
    fl.sort()
    idx = 0
    hist_idx=0
    r0.close()
    for fi in numpy.arange(len(fl)):
        f = fl[fi]
        res = h5py.File(f,'r')
        CS = CS + res["spec"][cc0,:,:]*numpy.conjugate(res["spec"][cc1,:,:])

        if (fi+1) % n_incoh_int == 0:

            if not quiet:
                print "%d/%d"%(idx,numpy.floor(len(fl)/n_incoh_int))
            plt.clf()
            #plt.subplot(211)
            rr = res["range"].value[0:200]/1e3

            fig, ax = plt.subplots()
            cax = ax.imshow(phase2imgcolors(CS[numpy.arange(200)[::-1],:]),aspect="auto",extent=[numpy.min(res["freq"]),numpy.max(res["freq"]), numpy.min(rr), numpy.max(rr)],cmap="hsv")
            plt.ylabel("Range (km)")
            plt.xlabel("Doppler (Hz)")
            plt.title("FM 92.3 MHz %s UTC\nRelative phase-Doppler-intensity image"%(stuffr.unix2datestr(res["t"].value)))
            norm = matplotlib.colors.Normalize(vmin=-numpy.pi, vmax=numpy.pi)

            cbar = fig.colorbar(cax, ticks=[0.0, 0.5, 1.0], cmap="hsv")
            cbar.ax.set_yticklabels(['$-\pi$', '0', '$\pi$'])
#            cb1 = matplotlib.colorbar.ColorbarBase(cax,
#                                                   norm=norm,
#                                                   orientation='vertical')

           # plt.subplot(212)
           # plt.pcolormesh(stuffr.comprz_dB(normalize(numpy.abs(CS[0:200,:])**2.0),fr=0.1))
           # plt.xlim([0,512])
           # plt.colorbar()
            plt.savefig("%s/%s/passive-%06d.png"%(dirn,pprefix,idx),orientation="landscape",dpi=150)
            plt.clf()
            idx = idx + 1
            CS[:,:]=0.0
        res.close()


# plot all individual spectra in directory
def analyze_all_passive_baselines(dirn,
                                  prefix="passive",
                                  pprefix="passive/plot",
                                  n_incoh_int=10,
                                  quiet=False):
    fl = glob.glob("%s/%s/cohin*.hdf5"%(dirn,prefix))
    r0 = h5py.File(fl[0],'r')
    n_channels = r0["spec"].shape[0]
    cs_idx = get_combinations(n_channels)
    CS = numpy.zeros( [cs_idx.shape[0],r0["spec"].shape[1],r0["spec"].shape[2]], dtype=numpy.complex64 )
    print(CS.shape)
    CS[:,:,:]=0.0
    SH = []
    os.system("mkdir -p %s/%s"%(dirn,pprefix))
    fl.sort()
    idx = 0
    hist_idx=0
    r0.close()
    for fi in numpy.arange(len(fl)):
        f = fl[fi]
        res = h5py.File(f,'r')
        CS = CS + res["spec"][cs_idx[:,0],:,:]*numpy.conjugate(res["spec"][cs_idx[:,1],:,:])

        if (fi+1) % n_incoh_int == 0:
            if not quiet:
                print "%d/%d"%(idx,numpy.floor(len(fl)/n_incoh_int))
            out_res = h5py.File("%s/%s/cs_spec-%06d.hdf5"%(dirn,pprefix,idx),'w')
            out_res["spec"]=CS
            out_res["range"]=res["range"]
            out_res["freq"]=res["freq"]
            idx = idx + 1
            CS[:,:,:]=0.0
            out_res.close()
            
        res.close()

# plot all individual spectra in directory
def plot_all_passive_gc(dirn,prefix="passive",pprefix="passive/plot",incoh_int=10):
    fl = glob.glob("%s/%s/cohin*.hdf5"%(dirn,prefix))
    r0 = h5py.File(fl[0],'r')
    print r0["clutter"].value.shape
    Nranges = len(r0["clutter"].value[0])
    Ntimes = math.floor(len(fl)/incoh_int) - 1
    R = numpy.zeros([Ntimes,Nranges],dtype=numpy.complex64)
    os.system("mkdir -p %s/%s"%(dirn,pprefix))
    fl.sort()
    idx = 0
    t0 = r0["t"].value
    r0.close()

    for fi in numpy.arange(len(fl)):
        f = fl[fi]
        res = h5py.File(f,'r')
        R[idx,:] = R[idx,:] + res["clutter"].value[0]*numpy.conjugate(res["clutter"].value[1])
        if (fi+1)%incoh_int == 0:
            print "%d/%d"%(idx,Ntimes)
            idx = idx + 1
        t1 = res["t"].value
        res.close()
        if idx+1 > Ntimes:
            break
    print "Done"
#    plt.pcolormesh(numpy.angle(numpy.transpose(R[:,numpy.arange(Nranges)[::-1]])))
#    plt.colorbar()

    fig, ax = plt.subplots()
    cax = ax.imshow(phase2imgcolors(numpy.transpose(R[:,numpy.arange(Nranges)[::-1]])),extent=[0,t1-t0, 0, Nranges*1.5],aspect="auto",cmap="hsv")
    plt.ylabel("Range (km)")
    plt.xlabel("Time (s)")
    plt.title("Ground clutter range time intensity relative phase difference")
    norm = matplotlib.colors.Normalize(vmin=-numpy.pi, vmax=numpy.pi)
    cbar = fig.colorbar(cax, ticks=[0.0, 0.5, 1.0], cmap="hsv")
    cbar.ax.set_yticklabels(['$-\pi$', '0', '$\pi$'])
    plt.savefig("%s/%s/gc.png"%(dirn,pprefix),orientation="landscape",dpi=150)
    plt.clf()
  #  res.close()

if __name__ == '__main__':

    parser = optparse.OptionParser()
    parser.add_option('-d', '--data_dir', dest='data_dir',
                      type='string', action='store', help='Full path to data directory')
    parser.add_option('-p', '--prefix', dest='prefix', type='string', action='store',
                      default="passive", help='Result directory prefix (default passive)')
    parser.add_option('-i', '--nint', dest='n_incoh_int', type='int', action='store',
                      default=10, help='Incoherent integration length (default 10)')
    parser.add_option('-0', '--channel_0', dest='channel_0', type='int', action='store', default=0,
                      help='Cross-correlation channel 0 (default %default)')
    parser.add_option('-1', '--channel_1', dest='channel_1', type='int', action='store', default=1,
                      help='Cross-correlation channel 1 (default %default)')
    parser.add_option('-q', '--quiet', dest='quiet', action='store_true',
                      help='No output')
    parser.add_option('-g', '--gc', dest='gc', action='store_true',
                      help='Plot ground clutter')

    (op, args) = parser.parse_args()

    if op.gc:
        plot_all_passive_gc(op.data_dir,
                            prefix="%s"%(op.prefix),pprefix="%s/plot"%(op.prefix))
    else:
        plot_all_passive(op.data_dir,
                         prefix="%s"%(op.prefix),pprefix="%s/plot"%(op.prefix),
                         n_incoh_int=op.n_incoh_int,dB_max=30.0,quiet=op.quiet,cc0=op.channel_0,cc1=op.channel_1)
