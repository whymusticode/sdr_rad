#!/usr/bin/env python
#
# Passive radar analysis
# incoherent integration using LPI
# (c) 2013, Juha Vierinen, MIT Haystack Observatory
#
import digital_rf_hdf5 as drf
import scipy, scipy.constants, scipy.linalg, numpy
import h5py
import stuffr
import os, time, optparse, errno, shutil, glob
import matplotlib.pyplot as plt
from multiprocessing import Process

import pyfftw

# create the theory matrix for a convolution (returns dictionary)
def create_theory_matrix(data,start_indx,max_range=200,coh_len=1000,tx_chan=None,rx_chan=None,idx_matrix=None):
    A = numpy.zeros([coh_len,max_range],dtype=numpy.complex64)
    z_tx = data.read_vector_c81d(start_indx-max_range,coh_len+max_range,tx_chan)
    z_rx = data.read_vector_c81d(start_indx-max_range,coh_len+max_range,rx_chan)
    A = z_tx[idx_matrix] # idx_matrix are indices gathering matrix from vector; runs fast!
    res = {"A":A,"m": z_rx[max_range + numpy.arange(coh_len)],"tx": z_tx[max_range + numpy.arange(coh_len)]}
    return(res)

# create index matrix for quick read of data into the theory matrix from a vector
def create_theory_idx_matrix(max_range=200,coh_len=1000):
    idx_matrix = numpy.zeros([coh_len,max_range],dtype=numpy.int)
    # build the theory matrix for convolution
    for i in range(coh_len):
        for j in range(max_range):
            idx_matrix[i,j] = i+max_range-j
    return(idx_matrix)

# ground clutter estimation using linear least squares (long coherence assumption); returns dictionary
def estimate_clutter(data,start_indx,coh_len=400,max_range=100,idx_matrix=None, tx_chan=0, rx_chan=1):
    r = create_theory_matrix(data,start_indx,max_range=max_range,coh_len=coh_len,idx_matrix=idx_matrix,
                             tx_chan=tx_chan, rx_chan=rx_chan)
    xhat = scipy.linalg.lstsq(r["A"],r["m"])[0]
    r["xhat"] = xhat
    return(r)

# estimate ground clutter, remove it, analyze spectrum using matched filtering (inverse filtering is ill-posed).
def analyze_passive_lpi(data,start_indx,coh_len=200,max_range=800,
                        idx_matrix_clutter=None, tx_chan=0, rx_chan=1,
                        ground_clutter_extent=70, lags=numpy.arange(50)+1, lag_decimation=10, n_incoh_samples=None):
    # initialize array for resulting echoes
    n_lags = len(lags)
    max_lag = lag_decimation*(1+numpy.max(lags))+1
    result_matrix = numpy.zeros([max_range,n_lags],dtype=numpy.complex64)
    # estimate ground clutter
    clutter = estimate_clutter(data,start_indx,
                               coh_len=coh_len+max_range+max_lag,
                               max_range=ground_clutter_extent,
                               tx_chan=tx_chan, rx_chan=rx_chan,
                               idx_matrix=idx_matrix_clutter)
    # remove clutter contribution from echoes
    m_without_clutter = clutter["m"][0:(coh_len+max_lag)]-numpy.dot(clutter["A"],clutter["xhat"])[0:(coh_len+max_lag)]
    tx = clutter["tx"]
    
    # sum freq domain toeplitz matrices for lpi
    tx_m = numpy.zeros([op.max_range,n_lags],dtype=numpy.complex64)
    tx_cov = numpy.zeros([op.max_range,n_lags],dtype=numpy.complex64)
#    window = stuffr.hanning(op.max_range)
    window = 1.0
    
    m_idx = numpy.arange(coh_len)
    m_without_clutter = clutter["m"][0:(coh_len+max_lag)]-numpy.dot(clutter["A"],clutter["xhat"])[0:(coh_len+max_lag)]
    for li in range(n_lags):
        for di in range(lag_decimation):
            lag_idx = di + lag_decimation*lags[li]+1
            mlp = m_without_clutter[m_idx]*numpy.conj(m_without_clutter[m_idx+lag_idx])
            tlp = tx[m_idx]*numpy.conj(tx[m_idx+lag_idx])
            for i in range(n_incoh_samples):
#                m0 = m_without_clutter[i*op.max_range + numpy.arange(max_range)]
#                m1 = m_without_clutter[i*op.max_range + numpy.arange(max_range)+lag_idx]
#                m_lp = numpy.fft.fft(m0*numpy.conj(m1))

    # initialize fftw
#    op.mlagp = pyfftw.n_byte_align_empty(op.max_range, 16, dtype='complex64')
#    op.mlagp[:]=0.0
#    op.txlagp = pyfftw.n_byte_align_empty(op.max_range, 16, dtype='complex64')
#    op.txlagp[:]=0.0
#    op.lpsum = pyfftw.n_byte_align_empty(op.max_range, 16, dtype='complex64')
#    op.lpsum[:]=0.0
#    op.fftmlagp = pyfftw.builders.fft(op.mlagp,threads=1)
#    op.ffttxlagp = pyfftw.builders.fft(op.txlagp,threads=1)
#    op.ifftlp = pyfftw.builders.ifft(op.lpsum,threads=1)


                op.mlagp[:] = mlp[i*op.max_range + numpy.arange(max_range)]
##                m_lp = numpy.fft.fft(mlp[i*op.max_range + numpy.arange(max_range)])
                m_lp = op.fftmlagp()
#                tx_lp = numpy.fft.fft(tx[i*op.max_range + numpy.arange(max_range)]*numpy.conj(tx[i*op.max_range + numpy.arange(max_range)+lag_idx]))
                op.txlagp[:]= tlp[i*op.max_range + numpy.arange(max_range)]
                tx_lp = op.ffttxlagp()
#                tx_lp = numpy.fft.fft(tlp[i*op.max_range + numpy.arange(max_range)])
                tx_m[:,li] = tx_m[:,li] + numpy.conj(tx_lp)*m_lp
                tx_cov[:,li] = tx_cov[:,li] + tx_lp*numpy.conj(tx_lp)
    for li in range(n_lags):
        result_matrix[:,li] = numpy.fft.ifft((1.0/tx_cov[:,li])*tx_m[:,li])
#    result_matrix = numpy.fft.ifft(tx_m)
#    stuffr.plot_cts(result_matrix)
    res = {"lag_matrix":result_matrix, "lags":lags*lag_decimation+1, "clutter": clutter["xhat"]}
    return(res)

# write an HDF5 file from a dictionary
def dict2hdf5(d,file_name):
    f = h5py.File(file_name,'w')
    for k in d.keys():
        f[k] = d[k]
    f.close()

# mkdir -p sucks with python (this was borrowed from stackoverflow.com)
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def get_lpi_results(dirn,dec=10):
    fl = glob.glob("%s/lpi-*.hdf5"%(dirn))
    fl.sort()
    n_files = len(fl)
    r0 = h5py.File(fl[0],"r")
    lagp = numpy.copy(r0["lag_matrix"])
    acf = numpy.zeros([lagp.shape[0],lagp.shape[1]*2],dtype=numpy.complex64)
    spec = numpy.zeros([lagp.shape[0],lagp.shape[1]*2],dtype=numpy.complex64)
    speca = numpy.zeros([lagp.shape[0],lagp.shape[1]*2],dtype=numpy.complex64)
    lagp[:,:]=0.0
    di = 0
    r0.close()
    var_est = 0.0
    scale_factor=0.0
    for i in range(n_files):
        print(i)
        r0 = h5py.File(fl[i],"r")
        lagm=numpy.copy(r0["lag_matrix"])
        print("spec")
        var_est = numpy.mean(numpy.abs(lagm[200:lagm.shape[0],:])**2.0)
        scale_factor = scale_factor + (1.0/var_est)
        lagp = lagp+(1.0/var_est)*lagm
        r0.close()
        if i % dec == 0 and i>0:
            lagp=lagp/scale_factor
            for ri in range(acf.shape[0]):
                acf[ri,0:lagp.shape[1]]=lagp[ri,:]
                acf[ri,(lagp.shape[1]):(2*lagp.shape[1])]=numpy.conj(lagp[ri,::-1])
                spec[ri,:]=numpy.fft.fftshift(numpy.fft.fft(acf[ri,:]))

            spec = spec/numpy.median(numpy.abs(spec))
            dB = stuffr.comprz_dB(spec[:,:])
            plt.pcolormesh(dB,vmin=0.0,vmax=30.0)
            plt.colorbar()
            plt.savefig("lpi_pr-%06d.png"%(i/dec),width=1024,height=768)
#            plt.show()
            plt.clf()
            lagp[:,:]=0.0
            scale_factor=0.0
            var_est=0.0

def get_lpi_results2(dirn,dec=10):
    fl = glob.glob("%s/lpi-*.hdf5"%(dirn))
    fl.sort()
    n_files = len(fl)
    r0 = h5py.File(fl[0],"r")
    lagp = numpy.copy(r0["lag_matrix"])
    acf = numpy.zeros([lagp.shape[0],lagp.shape[1]*2],dtype=numpy.complex64)
    spec = numpy.zeros([lagp.shape[0],lagp.shape[1]*2],dtype=numpy.complex64)
    speca = numpy.zeros([lagp.shape[0],lagp.shape[1]*2],dtype=numpy.complex64)
    lagp[:,:]=0.0
    di = 0
    r0.close()
    var_est = 0.0
    scale_factor=0.0
    for i in range(n_files):
        print(i)
        r0 = h5py.File(fl[i],"r")
        lagp=numpy.copy(r0["lag_matrix"])
        for ri in range(acf.shape[0]):
            acf[ri,0:lagp.shape[1]]=lagp[ri,:]
            acf[ri,(lagp.shape[1]):(2*lagp.shape[1])]=numpy.conj(lagp[ri,::-1])
            spec[ri,:]=numpy.fft.fftshift(numpy.fft.fft(acf[ri,:]))
        print("spec")
        var_est = numpy.mean(numpy.abs(spec[200:spec.shape[0],:])**2.0)
        scale_factor = scale_factor + (1.0/var_est)
        speca = speca+(1.0/var_est)*spec
        r0.close()
        if i % dec == 0 and i>0:
            speca = (1.0/scale_factor)*speca
            speca = speca/numpy.median(numpy.abs(speca))
#            dB = stuffr.comprz_dB(speca[0:200,:])
            dB = stuffr.comprz_dB(speca[:,:])
#            for si in range(spec.shape[1]):
 #               spec[:,si] = numpy.abs(spec[:,si])/numpy.mean(numpy.abs(spec[:,si]))
            plt.pcolormesh(dB,vmin=0.0,vmax=30.0)
            plt.colorbar()
            plt.savefig("lpi_pr-%06d.png"%(i/dec),width=1024,height=768)
#            plt.show()
            plt.clf()

            speca[:,:]=0.0
            scale_factor=0.0
            var_est=0.0

def plot_lpi_spec(f,dec=10):
    r = h5py.File(f,"r")
    spec = numpy.copy(r["lag_matrix"])
    for i in range(spec.shape[0]):
        spec[i,:]=10*numpy.log10(numpy.abs(numpy.fft.fftshift(numpy.fft.fft(r["lag_matrix"][i,:]))))
    r.close()
    return(spec)


        

def analyze_passive_radar(options):

    # start main processing

    op.data = drf.read_hdf5(op.data_dir)
    op.sample_rate = op.data.get_metadata(op.tx_chan)['sample_rate'].value

    op.result_dir = os.path.join(op.data_dir,op.rx_chan,op.result_dir_prefix)

    if op.reanalyze:
        if os.path.exists(op.result_dir):
            shutil.rmtree(op.result_dir)
    if not os.path.exists(op.result_dir):
        os.makedirs(op.result_dir)

    op.ranges = scipy.constants.c*numpy.arange(op.max_range)/op.sample_rate/2.0

    b = op.data.get_bounds(op.tx_chan)

    op.initial_sample_offset = op.max_range+1 + b[0]

    op.num_files = numpy.floor((b[1]-b[0])/(op.gc_coh_len))
    op.idx_matrix_clutter = create_theory_idx_matrix(max_range=op.ground_clutter_extent,coh_len=op.gc_coh_len+op.max_range+1+op.lag_decimation*(1+numpy.max(op.lags)))
    op.n_incoh_samples=int(op.gc_coh_len/op.max_range)
    op.n_existing_files = 0
    if not op.reanalyze:
        op.n_existing_files = len(glob.glob("%s/lpi-*.hdf5"%(op.result_dir)))

    # initialize fftw
    op.mlagp = pyfftw.n_byte_align_empty(op.max_range, 16, dtype='complex64')
    op.mlagp[:]=0.0
    op.txlagp = pyfftw.n_byte_align_empty(op.max_range, 16, dtype='complex64')
    op.txlagp[:]=0.0
    op.lpsum = pyfftw.n_byte_align_empty(op.max_range, 16, dtype='complex64')
    op.lpsum[:]=0.0
    op.fftmlagp = pyfftw.builders.fft(op.mlagp,threads=1)
    op.ffttxlagp = pyfftw.builders.fft(op.txlagp,threads=1)
    op.ifftlp = pyfftw.builders.ifft(op.lpsum,threads=1)
        
    compute_thread_n(0,op)
    # for ti in range(op.nthreads):
    #     p = Process(target=compute_thread_n,args=(ti,op,))
    #     p.start()
    #     procs.append(p)
    # for ti in range(op.nthreads):
    #     procs[ti].join()
    

# iterate through our own work units.
def compute_thread_n(n,op):
    print("compute %d"%(n))
    for i in numpy.arange(op.n_existing_files+n,op.num_files,op.nthreads*op.step):
        print("proc %d %d/%d"%(n,i,op.num_files))
        tic = time.time()
        t = op.data.get_metadata(op.tx_chan)["t0"].value+float(op.gc_coh_len*i)/op.sample_rate
        dname = stuffr.sec2dirname(t)
        fname = "%s/%s/lpi@%d.h5"%(op.result_dir,dname,t)
        if os.path.exists(fname):
            print("already exists, skipping")
        else:
            res = analyze_passive_lpi(op.data,
                                      start_indx=op.initial_sample_offset+op.gc_coh_len*i,
                                      coh_len=op.gc_coh_len,
                                      idx_matrix_clutter=op.idx_matrix_clutter,
                                      ground_clutter_extent=op.ground_clutter_extent,
                                      max_range=op.max_range,
                                      tx_chan=op.tx_chan, rx_chan=op.rx_chan,
                                      lags=op.lags, lag_decimation=op.lag_decimation, 
                                      n_incoh_samples=op.gc_coh_len/op.max_range)
            res["t"] = t
            res["range"] = op.ranges
        
            # pack output file
            os.system("mkdir -p %s/%s"%(op.result_dir,dname))
            dict2hdf5(res,"%s/%s/lpi@%d.h5"%(op.result_dir,dname,t))
            toc = time.time()
            if not op.quiet:
                print " File %d/%d dt %1.2f" % (i,op.num_files,toc-tic)


if __name__ == '__main__':

    parser = optparse.OptionParser()

    parser.add_option('-d', '--data_dir', dest='data_dir',
                      type='string', action='store',
                      help='Full path to data directory (no default)')
    parser.add_option('-c', '--gc_coh_len', dest='gc_coh_len', type='int', action='store',
                      default=102400, help='Ground clutter coherence length and number of samples to integrate (default %default)')
    parser.add_option('-l', '--lags', dest='lags', type='string', action='store',
                      default="numpy.arange(10)", help='Lag (default %default)')
    parser.add_option('-f', '--lag_decimation', dest='lag_decimation', type='int', action='store',
                      default=10, help='Lag decimation (default %default)')
    parser.add_option('-m', '--max_range', dest='max_range', type='int', action='store',
                      default=1024, help='Maximum range in samples (default %default)')
    parser.add_option('-g', '--ground_clutter_extent', dest='ground_clutter_extent',
                      type='int', action='store', default=70,
                      help='Ground clutter extent, in samples (default %default)')
    parser.add_option('-s', '--step', dest='step',
                      type='int', action='store', default=1,
                      help='Time step (default %default x)')
    parser.add_option('-i', '--initial_sample_offset', dest='initial_sample_offset',
                      type='int', action='store', default=0,
                      help='Initial offset into each data file, in samples (default %default)')
    parser.add_option('-p', '--result_dir_prefix', dest='result_dir_prefix',
                      type='string', action='store', default='passivew',
                      help='Results saved into directory with this name under data directory (default %default)')
    parser.add_option('-q', '--quiet', dest='quiet',
                      action="store_true", help='Quiet mode')
    parser.add_option('-x', '--reanalyze', dest='reanalyze',
                      action="store_true", help='Reanalyze')
    parser.add_option('-t', '--txchan', dest='tx_chan', type='string', action='store',
                      default="south", help='Transmitter channel (default %default)')
    parser.add_option('-r', '--rxchan', dest='rx_chan', type='string', action='store',
                      default="north", help='Receiver channel (default %default)')
    parser.add_option('-z', '--nthreads', dest='nthreads', type='int', action='store',
                      default=4, help='Number of threads (default %default)')
    (op, args) = parser.parse_args()
    op.lags=eval(op.lags)
    analyze_passive_radar(op)



