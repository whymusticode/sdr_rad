#!/usr/bin/env python
#
# Passive radar analysis
# (c) 2013, Juha Vierinen, MIT Haystack Observatory
#
import gdf
import scipy, scipy.constants, scipy.linalg, numpy
import h5py
import os, time, optparse, errno, shutil, glob
from multiprocessing import Process

# create the theory matrix for a convolution (returns dictionary)
def create_theory_matrix(data,start_indx,max_range=200,coh_len=1000,tx_chan=0,rx_chan=1,idx_matrix=None):
    A = numpy.zeros([coh_len,max_range],dtype=numpy.complex64)
    z_tx = gdf.read_vec(data,start_indx-max_range,coh_len+max_range,tx_chan)
    z_rx = gdf.read_vec(data,start_indx-max_range,coh_len+max_range,rx_chan)
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
def analyze_passive_mlfilt(data,start_indx,fft_len=512,coh_len=200,max_range=800,
                           idx_matrix_clutter=None,idx_matrix=None, tx_chan=0, rx_chan=1,
                           ground_clutter_extent=70):
    # initialize array for resulting echoes
    result_matrix = numpy.zeros([max_range,fft_len],dtype=numpy.complex64)
    # estimate ground clutter
    print("clutter %d"%(fft_len*coh_len+max_range))
    clutter = estimate_clutter(data,start_indx,
                               coh_len=fft_len*coh_len+max_range,
                               max_range=ground_clutter_extent,
                               tx_chan=tx_chan, rx_chan=rx_chan,
                               idx_matrix=idx_matrix_clutter)
    # remove clutter contribution from echoes
    m_without_clutter = clutter["m"][0:(fft_len*coh_len)]-numpy.dot(clutter["A"],clutter["xhat"])[0:(fft_len*coh_len)]
    # go through all time steps
    for i in range(fft_len):
        r = create_theory_matrix(data,start_indx+coh_len*i,max_range=max_range,coh_len=coh_len,idx_matrix=idx_matrix,
                                 tx_chan=tx_chan, rx_chan=rx_chan)
        A = r["A"]
        m = m_without_clutter[i*coh_len + numpy.arange(coh_len)]
        A_transpose = numpy.conjugate(numpy.transpose(A))
        # we have to do matched filtering here as a form of regularization because the
        # coherence time is too short for us to do maximum likelihood estimation
        # (not enough measurements). luckily the ground clutter is more coherent and
        # we can remove it's contribution to the sidelobes (see start of function)
        result_matrix[:,i] = numpy.dot(A_transpose,m)
    # estimate spectrum using simple boxcar window
    S = numpy.zeros([max_range,fft_len],dtype=numpy.complex64)
    for i in range(max_range):
        S[i,:] = numpy.fft.fftshift(numpy.fft.fft(result_matrix[i,:]))/fft_len
    res = {"spec":S, "clutter": clutter["xhat"]}
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

######################

if __name__ == '__main__':

    parser = optparse.OptionParser()
    parser.add_option('-c', '--coh_len', dest='coh_len', type='int', action='store',
                      default=200, help='Coherence length in samples (default %default)')
    parser.add_option('-n', '--fft_len', dest='fft_len', type='int', action='store',
                      default=512, help='FFT length (default %default)')
    parser.add_option('-m', '--max_range', dest='max_range', type='int', action='store',
                      default=800, help='Maximum range in samples (default %default)')
    parser.add_option('-g', '--ground_clutter_extent', dest='ground_clutter_extent',
                      type='int', action='store', default=70,
                      help='Ground clutter extent, in samples (default %default)')
    parser.add_option('-i', '--initial_sample_offset', dest='initial_sample_offset',
                      type='int', action='store', default=0,
                      help='Initial offset into each data file, in samples (default %default)')
    parser.add_option('-p', '--result_dir_prefix', dest='result_dir_prefix',
                      type='string', action='store', default='passive',
                      help='Results saved into directory with this name under data directory (default %default)')
    parser.add_option('-d', '--data_dir', dest='data_dir',
                      type='string', action='store',
                      help='Full path to data directory (no default)')
    parser.add_option('-q', '--quiet', dest='quiet',
                      action="store_true", help='Quiet mode')
    parser.add_option('-x', '--reanalyze', dest='reanalyze',
                      action="store_true", help='Reanalyze')
    parser.add_option('-t', '--txchan', dest='tx_chan', type='int', action='store',
                      default=0, help='Transmitter channel (default %default)')
    parser.add_option('-r', '--rxchan', dest='rx_chan', type='int', action='store',
                      default=1, help='Receiver channel (default %default)')
    parser.add_option('-z', '--nthreads', dest='nthreads', type='int', action='store',
                      default=1, help='Number of threads (default %default)')

    

    (op, args) = parser.parse_args()

    if op.initial_sample_offset < op.max_range:
        op.initial_sample_offset = op.max_range+1

    # start main processing

    op.data = gdf.new_gdf(op.data_dir)
    op.sample_rate = op.data['sample_rate']

    op.result_dir = os.path.join(op.data_dir,op.result_dir_prefix)

    if op.reanalyze:
        if os.path.exists(op.result_dir):
            shutil.rmtree(op.result_dir)
    if not os.path.exists(op.result_dir):
        os.makedirs(op.result_dir)

    op.freqs = numpy.linspace(-op.sample_rate/op.coh_len/2.0, op.sample_rate/op.coh_len/2.0, num=op.fft_len)

    op.ranges = scipy.constants.c*numpy.arange(op.max_range)/op.sample_rate/2.0
    op.num_files = numpy.round(op.data["max_n"]/(op.coh_len*op.fft_len))
    op.idx_matrix_clutter = create_theory_idx_matrix(max_range=op.ground_clutter_extent,coh_len=op.fft_len*op.coh_len+op.max_range)
    op.idx_matrix = create_theory_idx_matrix(max_range=op.max_range,coh_len=op.coh_len)
    op.n_existing_files = 0
    if not op.reanalyze:
        op.n_existing_files = len(glob.glob("%s/*.hdf5"%(op.result_dir)))

    for i in numpy.arange(op.n_existing_files,op.num_files):
        tic = time.time()
        res = analyze_passive_mlfilt(op.data,
                                     start_indx=op.initial_sample_offset+op.fft_len*op.coh_len*i,
                                     fft_len=op.fft_len,coh_len=op.coh_len,
                                     idx_matrix_clutter=op.idx_matrix_clutter,idx_matrix=op.idx_matrix,
                                     ground_clutter_extent=op.ground_clutter_extent,
                                     max_range=op.max_range,
                                     tx_chan=op.tx_chan, rx_chan=op.rx_chan)
        res["t"] = op.data["t0"]+float(op.fft_len*op.coh_len*i)/op.sample_rate
        res["freq"] = op.freqs
        res["range"] = op.ranges

        # pack output file
        dict2hdf5(res,os.path.join(op.result_dir,"cohint-%06d.hdf5" % i))
        toc = time.time()
        if not op.quiet:
            print " File %d/%d dt %1.2f" % (i,op.num_files,toc-tic)




