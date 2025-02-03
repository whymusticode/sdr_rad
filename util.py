
import numpy as np
from scipy.fft import fft, ifft
from scipy.sparse.linalg import lsqr
import time
import plotly.graph_objects as go

# from multiprocessing import Process, Queue
from gnuradio import gr, soapy, blocks
import numpy as np
import time
from numpy.lib.stride_tricks import as_strided


def quickPlot(thing):
    dbDat = np.log10(np.real(thing*np.conj(thing)))*10
    fig = go.Figure(data=go.Scatter(y=dbDat))
    fig.update_layout(template='plotly_dark',hovermode=False)
    fig.show()


def sdr_process_pluto(queue, frequency, sample_rate, uri=''):
    tb = gr.top_block()
    sdr = soapy.source('driver=plutosdr', "fc32", 1, uri, '', [''], [''])
    
    # Configure SDR
    sdr.set_sample_rate(0, sample_rate)
    sdr.set_frequency(0, frequency)
    sdr.set_gain(0, 20.0)
    
    # Create sink to capture data
    sink = blocks.vector_sink_c()
    tb.connect(sdr, sink)
    
    # Start flowgraph
    tb.start()
    
    while True:
        time.sleep(0.1)  # Collect for 100ms
        data = sink.data()
        sink.reset()
        if len(data) > 0:
            queue.put(data[:1024])  # Send 1024 samples

def sdr_process_RTL(queue, frequency, sample_rate, device_index=0):
    tb = gr.top_block()
    
    # Specify device by index
    dev_args = f'driver=rtlsdr,rtl={device_index}'
    
    sdr = soapy.source(dev_args, "fc32", 1, '',
                       'bufflen=16384', [''], [''])
    
    # Configure SDR
    sdr.set_sample_rate(0, sample_rate)
    sdr.set_frequency(0, frequency)
    sdr.set_gain(0, 'TUNER', 20.0)
    
    # Create sink to capture data
    sink = blocks.vector_sink_c()
    tb.connect(sdr, sink)
    
    # Start flowgraph
    tb.start()
    
    while True:
        time.sleep(0.1)
        data = sink.data()
        sink.reset()
        if len(data) > 0:
            queue.put(data[:1024])

def fftConv(x, taps, flag='valid'):
    N1 = len(x) - 1
    N2 = len(taps) - 1
    
    x = np.concatenate((x, np.zeros(N2, dtype=x.dtype)))
    taps = np.concatenate((taps, np.zeros(N1, dtype=x.dtype)))

    out = ifft(fft(x)*np.conj(fft(taps)))

    if flag == 'full':
        return np.roll(out, N2)
    if flag == 'same':
        return np.roll(out[:-N2], N2//2)
    elif flag == 'valid':
        return out[:-2*N2]
    else:
        raise Exception('undefined flag')

def db(x):
    dbx = np.log10(np.real(fft(x)*np.conj(fft(x))))*10
    return dbx - np.median(dbx)



def convMat(x, N):
    L = len(x) - N + 1
    strided = as_strided(x, shape=(L, N), 
                        strides=(x.strides[0], x.strides[0]),
                        writeable=False)
    return strided


def ax(x, tx, Ntaps, opt):
    N = len(tx)
    if opt == 'notransp':
        Nfft1 = N + Ntaps - 1
        y = ifft(fft(tx, Nfft1) * fft(x, Nfft1))
        return y[Ntaps-1:len(y)-Ntaps+1]
    else:
        Nfft2 = N + N + Ntaps - 1 - 1
        pre2 = fft(np.flipud(np.conj(tx)), Nfft2)
        y = ifft(pre2 * fft(x, Nfft2))
        return y[N-Ntaps:N]


def randn_iq(sz):
    return np.random.randn(sz).astype(np.complex64) + 1j*np.random.randn(sz).astype(np.complex64)


def db_atten(rxclean,rx): 
    return 10*np.log10(np.real(rxclean.conj().T@rxclean)/np.real(rx.conj().T@rx))



# N = 100000
# Ntaps = 100


# np.random.randn(Ntaps).astype(np.complex64)

# true_taps = randn_iq(Ntaps)
# tx = randn_iq(N) 
# A = convMat(tx, Ntaps)

# thing = A@np.conj(tx[100:200])
# quickPlot(thing)

# noise = randn_iq(N - Ntaps + 1)*1e-6
# rx = fftConv(tx,true_taps) + noise

# start = time.time()
# A = convMat(tx, Ntaps)
# make_mat_time = time.time() - start

# start = time.time()
# est_taps = np.conj(lsqr(A, rx,)[0]) # x0=true_taps 
# lineSolve_time = time.time() - start

# start = time.time()
# rx_est = fftConv(tx, est_taps, flag='valid')
# fft_time = time.time() - start

# start = time.time()
# A = convMat(tx, Ntaps)
# rx_est2 = A@np.conj(est_taps)
# mmult_time = time.time() - start

# rx_clean = rx - rx_est


# print(db_atten(rx_clean,rx))


# est_taps = randn_iq(Ntaps)
# for i in range(10):
#     rx_clean = rx - fftConv(tx,est_taps)
#     est_taps = est_taps + np.conj(A.conj().T@rx_clean*2e-7)  
#     print(db_atten(rx_clean,rx))

