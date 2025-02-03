#!/usr/bin/env python
#
# Passive radar analysis
# - multiple receiver channels
# (c) 2013, Juha Vierinen, MIT Haystack Observatory
#
import digital_rf_hdf5 as drf
import scipy, scipy.constants, scipy.linalg, numpy
import h5py
import os, time, optparse, errno, shutil, glob

show_progress = False

# create the theory matrix for a convolution (returns dictionary)
def create_theory_matrix(data,start_indx,max_range=200,coh_len=1000,tx_chan="tx",rx_chan=["rx"],idx_matrix=None):
    A = numpy.zeros([coh_len,max_range],dtype=numpy.complex64)
    z_tx = data.read_vector_c81d(start_indx-max_range,coh_len+max_range,tx_chan)
    z_rx = []
    for i in numpy.arange(len(rx_chan)):
        z_rx.append(data.read_vector_c81d(start_indx-max_range,coh_len+max_range,rx_chan[i])[max_range + numpy.arange(coh_len)])
    A = z_tx[idx_matrix] # idx_matrix are indices gathering matrix from vector; runs fast!
    res = {"A":A,"m": z_rx,"tx": z_tx[max_range + numpy.arange(coh_len)]}
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
def estimate_clutter(data,start_indx,coh_len=400,max_range=100,idx_matrix=None, tx_chan=0, rx_chan=[1]):
    r = create_theory_matrix(data,start_indx,max_range=max_range,coh_len=coh_len,idx_matrix=idx_matrix,
                             tx_chan=tx_chan, rx_chan=rx_chan)
    xhat = []
    for i in range(len(rx_chan)):
        xhat.append(scipy.linalg.lstsq(r["A"],r["m"][i])[0])
    r["xhat"] = xhat
    return(r)

# estimate ground clutter, remove it, analyze spectrum using matched filtering (inverse filtering is ill-posed).
def analyze_passive_mlfilt(data,start_indx,fft_len=512,coh_len=200,max_range=800,
                           idx_matrix_clutter=None,idx_matrix=None, tx_chan=0, rx_chan=[1],
                           ground_clutter_extent=70):
    # initialize array for resulting echoes
    n_rx_channels = len(rx_chan)
    result_matrix = numpy.zeros([n_rx_channels,max_range,fft_len],dtype=numpy.complex64)
    if show_progress:
        print("estimating clutter")
    # estimate ground clutter for all channels, and remove from measurements
    clutter = estimate_clutter(data,start_indx,
                               coh_len=fft_len*coh_len+max_range,
                               max_range=ground_clutter_extent,
                               tx_chan=tx_chan, rx_chan=rx_chan,
                               idx_matrix=idx_matrix_clutter)
    m_without_clutter = []
    for i in range(n_rx_channels):
        m_without_clutter.append(clutter["m"][i][0:(fft_len*coh_len)]-numpy.dot(clutter["A"],clutter["xhat"][i])[0:(fft_len*coh_len)])
    # go through all time steps
    if show_progress:
        print("doing matched filter")
    for i in range(fft_len):
        if show_progress:
            print("doing matched filter %d/%d"%(i,fft_len))
        r = create_theory_matrix(data,start_indx+coh_len*i,max_range=max_range,coh_len=coh_len,idx_matrix=idx_matrix,
                                 tx_chan=tx_chan, rx_chan=rx_chan)
        for j in range(n_rx_channels):
            A = r["A"]
            m = m_without_clutter[j][i*coh_len + numpy.arange(coh_len)]
            A_transpose = numpy.conjugate(numpy.transpose(A))
            # we have to do matched filtering here as a form of regularization because the
            # coherence time is too short for us to do maximum likelihood estimation
            # (not enough measurements). luckily the ground clutter is more coherent and
            # we can remove it's contribution to the sidelobes (see start of function)
            result_matrix[j,:,i] = numpy.dot(A_transpose,m)

    # estimate spectrum using simple boxcar window
    S = numpy.zeros([n_rx_channels,max_range,fft_len],dtype=numpy.complex64)
    for j in range(n_rx_channels):
        for i in range(max_range):
            S[j,i,:] = numpy.fft.fftshift(numpy.fft.fft(result_matrix[j,i,:]))/fft_len
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
                      default=500, help='FFT length (default %default)')

    parser.add_option('-m', '--max_range', dest='max_range', type='int', action='store',
                      default=800, help='Maximum range in samples (default %default)')

    parser.add_option('-g', '--ground_clutter_extent', dest='ground_clutter_extent',
                      type='int', action='store', default=70,
                      help='Ground clutter extent, in samples (default %default)')

    parser.add_option('-i', '--initial_sample_offset', dest='initial_sample_offset',
                      type='int', action='store', default=0,
                      help='Initial offset into each data file, in samples (default %default)')

    parser.add_option('-p', '--result_dir', dest='result_dir',
                      type='string', action='store', default='/data1/results/fm_south105.1/passive',
                      help='Results saved into directory with this name under data directory (default %default)')

    parser.add_option('-d', '--data_dir', dest='data_dir',
                      type='string', action='store',
                      help='Full path to data directory (no default)')

    parser.add_option('-q', '--quiet', dest='quiet',
                      action="store_true", help='Quiet mode')

    parser.add_option('-x', '--reanalyze', dest='reanalyze',
                      action="store_true", help='Reanalyze')

    parser.add_option('-t', '--txchan', dest='tx_chan', type='string', action='store',
                      default="ch000", help='Transmitter channel (default %default)')

    parser.add_option('-r', '--rxchan', dest='rx_chan', type='string', action='store',
                      default="[\"ch001\",\"ch002\",\"ch003\",\"ch004\",\"ch005\",\"ch006\",\"ch007\"]",
                      help='Receiver channel (default %default)')

    (op, args) = parser.parse_args()

    if op.initial_sample_offset < op.max_range:
        op.initial_sample_offset = op.max_range+1
    op.rx_chan = eval(op.rx_chan)

    # start main processing
    data = drf.read_hdf5(op.data_dir)
    channels = data.get_channels()
    metadata = data.get_metadata(channels[0])
    sample_rate = metadata['sample_rate'].value

    result_dir = op.result_dir

    if op.reanalyze:
        if os.path.exists(result_dir):
            shutil.rmtree(result_dir)
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    freqs = numpy.linspace(-sample_rate/op.coh_len/2.0, sample_rate/op.coh_len/2.0, num=op.fft_len)

    ranges = scipy.constants.c*numpy.arange(op.max_range)/sample_rate/2.0

    bounds = data.get_bounds(channels[0])
    num_files = numpy.round((bounds[1]-bounds[0])/(op.coh_len*op.fft_len))
    
    idx_matrix_clutter = create_theory_idx_matrix(max_range=op.ground_clutter_extent,coh_len=op.fft_len*op.coh_len+op.max_range)
    idx_matrix = create_theory_idx_matrix(max_range=op.max_range,coh_len=op.coh_len)
    
    # analyze whole directory
    for i in numpy.arange(0,num_files):
        if not op.quiet:
            print " File %d/%d" % (i,num_files)
        res_sec = int((i*op.fft_len*op.coh_len)/sample_rate) + bounds[0]/sample_rate

        if os.path.isfile(os.path.join(result_dir,"cohint-%d.hdf5" % (res_sec))):
            print("file exists")
        else:
            res = analyze_passive_mlfilt(data,
                                         start_indx=op.initial_sample_offset+op.fft_len*op.coh_len*i + bounds[0],
                                         fft_len=op.fft_len,coh_len=op.coh_len,
                                         idx_matrix_clutter=idx_matrix_clutter,idx_matrix=idx_matrix,
                                         ground_clutter_extent=op.ground_clutter_extent,
                                         max_range=op.max_range,
                                         tx_chan=op.tx_chan, rx_chan=op.rx_chan)
            
            res["t"] = bounds[0]/sample_rate+float(op.fft_len*op.coh_len*i)/sample_rate
            res["freq"] = freqs
            res["range"] = ranges

            # pack output file
            dict2hdf5(res,os.path.join(result_dir,"cohint-%d.hdf5" % (res_sec)))


