#!/usr/bin/python
#
# plotting routines for passive radar results
#
if __name__ == '__main__':
    import matplotlib
    matplotlib.use("Agg")

import stuffr
import digital_rf_hdf5 as drf
import glob, os, numpy, h5py
import matplotlib.pyplot as plt
import os, time, optparse, errno

# mkdir -p sucks with python (this was borrowed from stackoverflow.com)
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

# plot all individual spectra in directory
def get_gc(dirn,prefix="passive"):
    fl = glob.glob("%s/%s/cohin*.hdf5"%(dirn,prefix))
    fl.sort()
    n_files = (len(fl)-10)
    print(n_files)
    res = numpy.zeros([n_files],dtype=numpy.complex64)
    for i in range(n_files):
        r = h5py.File(fl[i],'r')
        res[i] = r["clutter"][0]
        if i % 100 == 0:
            print(".")
        r.close()
    return(res)

# plot all individual spectra in directory
def plot_all_passive_rti(dirn,prefix="passive",pprefix="passive/plot",incoh_int=10, dB_max=30.0,Ntimes=None):
    fl = glob.glob("%s/%s/cohin*.hdf5"%(dirn,prefix))
    r0 = h5py.File(fl[0],'r')
    Nranges = r0["spec"].shape[0]
    print("Nranges %d"%(Nranges))
    if Ntimes is None:
        Ntimes = int(numpy.floor(len(fl)/incoh_int) - 1)
    print("n_times %d"%(Ntimes))
    R = numpy.zeros([Ntimes,Nranges],dtype=numpy.complex64)
    S = numpy.array(numpy.copy(r0["spec"]),dtype=numpy.float32)
    S[:,:]=0.0
    os.system("mkdir -p %s/%s"%(dirn,pprefix))
    fl.sort()
    idx = 0
    t0 = r0["t"].value
    times = numpy.zeros([Ntimes])

    for fi in numpy.arange(len(fl)):
        f = fl[fi]
        res = h5py.File(f,'r')
        S = S + numpy.abs(res["spec"])**2.0
        
        if (fi+1)%incoh_int == 0:
            print "%d/%d"%(idx,Ntimes)
            times[idx]=res["t"].value
            for ri in range(Nranges):
                R[idx,ri] = numpy.max(S[ri,:])
            idx = idx + 1
            S[:,:]=0.0
        t1 = res["t"].value
        res.close()
        if idx+1 > Ntimes:
            break
    print "Done"
    for ti in range(Ntimes):
        R[ti,:]=R[ti,:]/numpy.median(R[ti,:])
    dBP = stuffr.comprz_dB(R[0:Ntimes,:])
    plt.clf()
    plt.pcolormesh(numpy.transpose(dBP),cmap="gist_yarg")
    plt.colorbar()
    plt.savefig("%s/%s/rti.jpg"%(dirn,pprefix),orientation="landscape",dpi=200)
    plt.clf()
    plt.pcolormesh(times, r0["range"].value/1e3,numpy.transpose(dBP),cmap="gist_yarg")
    plt.ylim([0.0,400.0])
    plt.xlim([numpy.min(times),numpy.max(times)])
    plt.colorbar()
    plt.savefig("%s/%s/rti_airplanes.jpg"%(dirn,pprefix),orientation="landscape",dpi=200)
    plt.clf()
    r0.close()

# plot all individual spectra in directory
def plot_all_passive(dirn,prefix="passive",pprefix="passive/plot",n_incoh_int=10,dB_max=30.0, quiet=False, replot=False, op=False):
    g = drf.read_hdf5(dirn)
    fuck_this = "%s/%s/cohin*.hdf5"%(dirn,prefix)
    print(fuck_this)
    fl = glob.glob(fuck_this)
    r0 = h5py.File(fl[0],'r')
    S = numpy.copy(numpy.abs(r0["spec"][0,:,:])**2.0)
    r0.close()
    S[:,:]=0.0
    SH = []
    os.system("mkdir -p %s/%s"%(dirn,pprefix))
    fl.sort()
    idx = 0
    hist_idx=0
    clutter_max=0.0
    n_plots = 0
    if not replot:
        n_plots = len(glob.glob("%s/%s/plot/*.jpg"%(dirn,prefix)))
        idx = n_plots

    for fi in numpy.arange(n_plots*n_incoh_int,len(fl)):
        f = fl[fi]
        res = h5py.File(f,'r')
        clutter_max = clutter_max+numpy.max(numpy.abs(res["clutter"].value)**2)
        S = S + numpy.abs(res["spec"][0,:,:])**2.0

        if (fi+1) % n_incoh_int == 0:
            if not quiet:
                print "%d/%d"%(idx,numpy.floor(len(fl)/n_incoh_int))
            plt.clf()
            dBP = 10.0*numpy.log10(numpy.abs(S))
            dmin = numpy.median(dBP)
            snr=numpy.max(dBP)-dmin
            cnr = 10.0*numpy.log10(clutter_max/numpy.median(numpy.abs(S)))
            print "CNR %1.2f SNR %1.2f"%(cnr,snr)
            dBP = dBP-dmin
            plt.pcolormesh(res["freq"].value, res["range"].value/1e3, dBP,vmin=0.0,vmax=dB_max,cmap="gist_yarg")
            plt.xlim([numpy.min(res["freq"]),numpy.max(res["freq"])])
#            plt.ylim([0.0,numpy.max(res["range"])/1e3])
            plt.ylim([0.0,numpy.max(res["range"].value/1e3)])
            plt.ylabel("Range (km)")
            plt.xlabel("Doppler (Hz)")
            plt.title("FM %1.2f MHz, Passive Radar, %s\nCNR=%01.2f dB, Max SNR=%01.2f dB"%(g.get_metadata(op.channel)["center_frequencies"].value/1e6,
                                                                                           stuffr.unix2datestr(res["t"].value),cnr,snr))
            plt.colorbar()

            plt.savefig("%s/%s/passive-%06d.jpg"%(dirn,pprefix,idx),orientation="landscape",dpi=200)
            plt.clf()
            idx = idx + 1
            S[:,:]=0.0
            clutter_max=0.0
        res.close()

if __name__ == '__main__':

    parser = optparse.OptionParser()
    parser.add_option('-d', '--data_dir', dest='data_dir',
                      type='string', action='store', help='Full path to data directory')
    parser.add_option('-p', '--prefix', dest='prefix', type='string', action='store',
                      default="passive", help='Result directory prefix (default passive)')
    parser.add_option('-i', '--nint', dest='n_incoh_int', type='int', action='store',
                      default=10, help='Incoherent integration length (default 10)')
    parser.add_option('-q', '--quiet', dest='quiet', action='store_true',
                      help='No output')
    parser.add_option('-r', '--replot', dest='replot', action='store_true',
                      help='Replot all images (default False).')
    parser.add_option('-c', '--channel', dest='channel', action='store',type="string",default="fm_south105.10",
                      help='Channel (default %default).')

    (op, args) = parser.parse_args()
#    plot_all_passive_rti(op.data_dir,
#                         prefix="%s"%(op.prefix),pprefix="%s/plot"%(op.prefix),Ntimes=None,
#                         incoh_int=op.n_incoh_int,dB_max=30.0)

    plot_all_passive(op.data_dir,
                     prefix="%s/%s"%(op.channel,op.prefix),pprefix="%s/%s/plot"%(op.channel,op.prefix),
                     n_incoh_int=op.n_incoh_int,dB_max=30.0,quiet=op.quiet, replot=op.replot, op=op)
