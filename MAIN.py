
# https://pysdr.org/content/pluto.html
# https://wiki.analog.com/university/tools/pluto
# https://www.youtube.com/watch?v=ph0Kv4SgSuI # dual receive X dual transmit

# XC7Z010-1CLG225C4334
# https://github.com/analogdevicesinc/hdl/tree/main/projects/pluto

# on mac, connecting is a little different (a little fuckign stupid): 
# ls -l /dev/tty.* 
# crw-rw-rw-  1 root  wheel   17,   2 Nov  7 15:28 /dev/tty.usbmodem1414
# screen /dev/tty.usbmodem1414 115200

# XC7Z010-1CLG225

# https://wiki.analog.com/university/tools/pluto/hacking/hardware

# fw_setenv compatible ad9364
# reboot

# iio_info -s

import adi
import time
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import tensorflow as tf
import util
import importlib
importlib.reload(util)
from numpy.lib.stride_tricks import as_strided
import plotly.graph_objects as go
import sledge


fft = np.fft.fft
ifft = np.fft.ifft
fftshift = np.fft.fftshift
conj = np.conj
log10 = np.log10
abs = np.abs
real = np.real
flipud = np.flipud
from scipy.linalg import lstsq


sdr = adi.ad9361(uri='ip:192.168.2.1')

# fs = 30.72e6 #  >= 521e3 and <=30.72e6 if 2 channels 
# fc = 98e6  # 104.1e6
# Tint = .01
# Nsamp = int(Tint*fs) #samples per buffer.  Can be different for Rx and Tx


fs = 521e3 #  >= 521e3 and <=30.72e6 if 2 channels 
fc = 89.2e6  # 104.1e6
# for fc in range(int(88e6),int(108e6),int(2e6)):
Tint = 2
Nsamp = int(Tint*fs) #samples per buffer.  Can be different for Rx and Tx

rx_mode = "manual"  # can be "manual" or "slow_attack"


sdr.rx_enabled_channels = [0, 1]
sdr.sample_rate = int(fs)
sdr.fc = int(fc)
sdr.gain_control_mode = rx_mode
sdr.rx_buffer_size = int(Nsamp)
sdr.rx_hardwaregain_chan0 = int(-10)
sdr.rx_hardwaregain_chan1 = int(-10)
sdr.rx_rf_bandwidth = int(200e3)

# rxdata = [x.astype(np.complex64) for x in sdr.rx()]
# freq = fftshift(np.fft.fftfreq(Nsamp, 1/fs)) + fc
# spectrum = np.fft.fftshift(np.fft.fft(rxdata[1] ))
# powerDb = 20*np.log10(np.abs(spectrum))
# fig = go.Figure()
# fig.add_trace(go.Scatter(x = freq,y = powerDb))
# fig.show()




RXidx = 0
TXidx = 1
N = 100
Ntaps = 30
stop = Ntaps//2
start = stop-1
p = {'fc':fc,'fs':fs,'vMax':300,'rMax':200e3}
frames = []
for r in range(2):  
    # tic = time.time()
    rxdata = [x.astype(np.complex64) for x in sdr.rx()]
    # toc = time.time() - tic
    # print(Tint/toc)

    maxOffset = 100000 
    spectrum1 = fft(rxdata[0][0:maxOffset]) 
    spectrum2 = fft(rxdata[1][0:maxOffset]) 
    xcorr = fftshift(ifft(spectrum1*conj(spectrum2))) 

    xcorr_dB = 20*log10(abs(xcorr)) 
    spec_dB = 20*log10(fftshift(abs(spectrum1))) 
    spec_dB2 = 20*log10(fftshift(abs(spectrum2))) 

    fig = make_subplots(rows=3, cols=1, subplot_titles=( 'signal','Spectrum','Cross-correlation')) 
    _=fig.add_trace(go.Scatter(y=real(rxdata[0][0:1000])), row=1, col=1) 
    _=fig.add_trace(go.Scatter(y=real(rxdata[1][0:1000])), row=1, col=1) 

    _=fig.add_trace(go.Scatter(y=spec_dB), row=2, col=1) 
    _=fig.add_trace(go.Scatter(y=spec_dB2), row=2, col=1) 
    _=fig.add_trace(go.Scatter(y=xcorr_dB), row=3, col=1) 
    _=fig.update_layout(title_text="Channel Coherence Analysis", height=800) 
    fig.show() 
    

    b = rxdata[RXidx][start:-stop]
    A_batched = tf.convert_to_tensor(as_strided(rxdata[TXidx], 
    shape=(len(b)//N, N, Ntaps), 
    strides=(rxdata[TXidx].strides[0]*N, rxdata[TXidx].strides[0], rxdata[TXidx].strides[0])))
    b_batched = tf.convert_to_tensor(b[:len(b)//N*N].reshape(-1,N,1))
    try:
        x = tf.linalg.lstsq(A_batched, b_batched)
    except ValueError as e:
        print(f"Error: {e}")

    cleaned_batched = b_batched - tf.matmul(A_batched, x)
    cleaned_flat = tf.reshape(cleaned_batched, [-1])  # Back to 1D
    rxClean = tf.concat([tf.zeros(start, dtype=cleaned_flat.dtype),cleaned_flat,tf.zeros(stop, dtype=cleaned_flat.dtype)], axis=0)
    rxClean = np.array(rxClean)
    print(fc)
    print(util.db_atten(rxClean,b))

    range_vals, velocity, image = sledge.sledge(rxdata[0], rxClean, p)
    z = 10*np.log10(image.real**2 + image.imag**2)

    fig = go.Figure(data=go.Heatmap(z=z,x=range_vals/1e3,y=velocity,))
    _ = fig.update_layout(template='plotly_dark',hovermode=False)
    fig.show()




    # frames.append(go.Frame(data=[go.Heatmap(z=z, x=range_vals/1e3, y=velocity)]))








# rxdata[0] = carrier
# demodulated, carrier = fm_demod(rxdata[1], fs, fc)
# rxdata[1] = carrier



# N = 100
# Ntaps = 10
# data = util.randn_iq(N)
# taps = util.randn_iq(Ntaps)
# A = util.convMat(data,Ntaps)
# # b = A@taps + util.randn_iq(N-Ntaps+1)*.0001
# b = util.fftConv(data,taps) + util.randn_iq(N-Ntaps+1)*.0001
# x = lstsq(A, b)[0]
# rxClean = np.concat((np.zeros(Ntaps//2-1), b - A@x,np.zeros(Ntaps//2)))
# util.db_atten(rxClean,b)


# asdf = 20*log10(abs(fftshift(fft(rxClean))))
# util.quickPlot(fftshift(fft(rxClean)))

# util.convMat(np.array(range(10)),5)





# import numpy as np
# from pyapril import RDTools, detector, channelPreparation

# # # Load your IQ data
# # iq_data = np.load('your_data.npy')
# ref_ch = rxdata[0]
# surv_ch = rxdata[1]

# # Apply channel preparation
# tic = time.time()
# surv_ch_filtered = channelPreparation.time_domain_filter_surveillance(
#     ref_ch, 
#     surv_ch, 
#     "ECA",  # Filter method
#     K=64,   # Filter dimension
#     D=2     # Decimation factor
# )
# toc = time.time()-tic

# # Generate range-Doppler matrix
# rd_matrix = detector.cc_detector_ons(
#     ref_ch,
#     surv_ch_filtered,
#     200e3,  # Sample rate
#     100,    # Max Doppler
#     128,    # Number of range bins
#     verbose=0
# )

# # Visualize results
# RDTools.export_rd_matrix_img(
#     'output.png',
#     rd_matrix,
#     max_Doppler=100,
#     dyn_range=40,
#     dpi=800
# )




# tic = time.time()
# rxClean = np.zeros(len(rxdata[0]), dtype=complex)
# for i in range(len(b)//N):
#     idx = slice(i*N, (i+1)*N)  # Fixed indexing
#     Atmp = A[idx,:]
#     btmp = b[idx]
#     x = np.linalg.lstsq(Atmp, btmp, rcond=None)[0]
#     rxClean[start + i*N:start + (i+1)*N] = btmp - Atmp@x
# rxClean = np.concatenate((np.zeros(start), rxClean[start:len(rxdata[0])-stop], np.zeros(stop)))
# toc = time.time()-tic
# util.db_atten(rxClean,b)

# import tensorflow as tf

# batched_A = tf.reshape(A, [-1, N, Ntaps]) 
# batched_b = tf.reshape(b, [-1, N])
# X = tf.linalg.lstsq(batched_A, tf.expand_dims(batched_b, -1))

# # Get cleaned signal by doing batch matrix multiply
# rxClean_batched = batched_b - tf.squeeze(tf.matmul(batched_A, tf.expand_dims(X, -1)))

# # Reshape back to original dimensions
# rxClean[start:len(rxdata[0])-stop] = tf.reshape(rxClean_batched, [-1])


# import tensorflow as tf
# A = tf.random.normal([100, 5, 5])  # 100 matrices, each 5x5 
# b = tf.random.normal([100, 5, 1])  # 100 vectors, each 5x1
# x = tf.linalg.solve(A, b)  # x shape: [100, 5, 1]




# import numpy as np
# import cupy as cp
# from concurrent.futures import ProcessPoolExecutor
# import time
# from functools import partial

# def process_batch_gpu(batch_data, A_gpu, start):
#     """Process a single batch on GPU"""
#     idx, btmp = batch_data
#     Atmp_gpu = cp.array(A[idx,:])
#     btmp_gpu = cp.array(btmp)
    
#     # Solve least squares on GPU
#     x = cp.linalg.lstsq(Atmp_gpu, btmp_gpu, rcond=None)[0]
    
#     # Compute cleaned signal
#     cleaned = btmp_gpu - Atmp_gpu @ x
    
#     return idx, cp.asnumpy(cleaned)

# def optimized_clean_signal(rxdata, Ntaps, N, n_workers=32, batch_size=1000):
#     """
#     Clean signal using GPU-accelerated parallel processing
    
#     Args:
#         rxdata: Input signal data
#         Ntaps: Number of taps
#         N: Processing window size
#         n_workers: Number of parallel workers
#         batch_size: Size of batches for parallel processing
#     """
#     tic = time.time()
    
#     # Initialize output array
#     rxClean = np.zeros(len(rxdata[0]), dtype=complex)
    
#     # Compute start/stop
#     start = Ntaps - 1
#     stop = Ntaps - 1
    
#     # Create convolution matrix
#     A = util.convMat(rxdata[1], Ntaps)
#     b = rxdata[0][start:-stop]
    
#     # Transfer constant data to GPU
#     A_gpu = cp.array(A)
    
#     # Create batches
#     total_batches = len(b) // N
#     batch_indices = [(slice(i*N, (i+1)*N), b[i*N:(i+1)*N]) 
#                     for i in range(total_batches)]
    
#     # Process in parallel across multiple GPUs
#     process_func = partial(process_batch_gpu, A_gpu=A_gpu, start=start)
    
#     with ProcessPoolExecutor(max_workers=n_workers) as executor:
#         # Process batches in parallel
#         for idx, cleaned_batch in executor.map(process_func, batch_indices):
#             rxClean[start + idx.start:start + idx.stop] = cleaned_batch
    
#     # Combine results
#     rxClean = np.concatenate((
#         np.zeros(start), 
#         rxClean[start:len(rxdata[0])-stop],
#         np.zeros(stop)
#     ))
    
#     toc = time.time() - tic
#     print(f"Processing time: {toc:.2f} seconds")
    
#     return rxClean

# # For multi-GPU support
# def get_gpu_id():
#     """Get GPU ID for worker"""
#     worker_id = multiprocessing.current_process()._identity[0]
#     return worker_id % torch.cuda.device_count()

# # Set batch size based on GPU memory
# BATCH_SIZE = 1000
# N_WORKERS = 32  # Adjust based on number of CPU cores

# cleaned_signal = optimized_clean_signal(
#     rxdata,
#     Ntaps=Ntaps,
#     N=N,
#     n_workers=N_WORKERS,
#     batch_size=BATCH_SIZE
# )














# import numpy as np
# from scipy import linalg
# from concurrent.futures import ProcessPoolExecutor

# def mm3d(a, b):
#     """3D matrix multiplication matching MATLAB's behavior
#     For matrices a(i,j,k) and b(l,m,k), returns c(i,m,k)
#     """
#     return np.einsum('ijk,jlk->ilk', a, b)

# def process_block(args):
#     """Process single block for parallel execution"""
#     H_block, Hp_block, RX_block = args
#     HH = mm3d(Hp_block, H_block)
#     tmp = np.linalg.solve(HH, Hp_block)
#     return mm3d(tmp, RX_block)

# # def remove_direct_vect(win_len, n_taps, tx, rx, n_workers=4):

# win_len = 100
# n_taps = 20
# tx = rxdata[0]
# rx = rxdata[1]

# n = len(tx)
# n_wins = n // win_len

# # Create FIR filter matrix
# input_sig = tx
# H_pre = np.tile(input_sig, (n_taps, 1))
# cols = (H_pre.shape[1]) // (n + 1)
# n_T = (cols - 1) // 2  # actual number of taps

# # Build H matrix
# HH = np.vstack([input_sig[-(n_T-1):], H_pre[:((n+1)*cols-n_T)]])
# HH = HH.reshape((n+1, cols))
# HH = HH[:-1]

# # Reshape to 3D
# H = HH.reshape((win_len, n_wins, cols))
# H = np.transpose(H, (0, 2, 1))

# # Conjugate transpose
# Hp = np.conj(np.transpose(H, (1, 0, 2)))

# # Reshape rx for processing
# RX2 = rx.reshape((win_len, 1, n_wins))

# # Process blocks in parallel
# blocks = [(H[:,:,i:i+1], Hp[:,:,i:i+1], RX2[:,:,i:i+1]) 
#             for i in range(n_wins)]

# with ProcessPoolExecutor(max_workers=n_workers) as executor:
#     ws_blocks = list(executor.map(process_block, blocks))

# # Combine results
# ws = np.concatenate(ws_blocks, axis=2)

# # Calculate direct path
# direct = mm3d(H, ws)

# # Remove direct path
# rx_clean = rx - direct.flatten()




# remove_direct_vect(100,20,rxdata[0],rxdata[1])


# a = np.random.random((10,5,5))  # (i,j,k)
# b = np.random.random((5,3,5))   # (l,j,k)
# c = mm3d(a, b)                  # result is (i,l,k) = (10,3,5)
# print(c.shape)  # Should print (10, 3, 5)


# a[:,:,0]@b[:,:,0].T
# c[:,:,0]

# winLen = 100;
# nTaps = 20;

# N = 1000
# Ntaps = 20
# stop = Ntaps//2
# start = stop-1
# A = util.convMat(rxdata[1][0:N],Ntaps)
# b = rxdata[0][0:N][start:-stop]
# x = lstsq(A, b)[0]
# rxClean = np.concat((np.zeros(start), b - A@x,np.zeros(stop)))
# util.db_atten(rxClean,b)




# np.concat(( np.array((0.,0.)),np.array((1.,1.)) ))

# def Quick_plot(data):
#     powerDb = 10 * log10(abs(data)**2)
#     fig = go.Figure()
#     Nds = int(round(len(powerDb)/1e3))
#     fig.add_trace(go.Scatter(x = freq[::Nds],y = powerDb[::Nds]))
#     fig.show()

# freq = np.fft.fftshift(np.fft.fftfreq(Nsamp, 1/fs)) + fc
# spectrum = np.fft.fftshift(np.fft.fft(rxdata[1] ))
# powerDb = 10 * np.log10(np.abs(spectrum)**2)

# fig = go.Figure()
# Nds = int(round(len(powerDb)/10e3))
# fig.add_trace(go.Scatter(x = freq[::Nds],y = powerDb[::Nds]))
# fig.show()
# print(f"Total time: {time.time() - tic} seconds")


####### don't fucking touch these 
# sdr.rx_lo = int(fc)
# sdr.tx_lo = int(fc)


###### no need to store data 
# specData = np.vstack([specData[1:], powerDb])
# toc = time.time() - tic
# print(Texpected/toc) # check that we aren't loosing that much data 


### DSP from the original example 
# Rx_0=data[0]
# Rx_1=data[1]
# Rx_total = Rx_0 + Rx_1
# NumSamples = len(Rx_total)
# win = np.hamming(NumSamples)
# y = Rx_total * win
# sp = np.absolute(np.fft.fft(y))
# sp = sp[1:-1]
# sp = np.fft.fftshift(sp)
# s_mag = np.abs(sp) / (np.sum(win)/2)    # Scale FFT by window and /2 since we are using half the FFT spectrum
# s_dbfs = 20*np.log10(s_mag/(2**12))     # Pluto is a 12 bit ADC, so use that to convert to dBFS
# xf = np.fft.fftfreq(NumSamples, ts)
# xf = np.fft.fftshift(xf[1:-1])/1e6
# plt.plot(xf, s_dbfs)
# plt.xlabel("frequency [MHz]")
# plt.ylabel("dBfs")
# plt.draw()
# plt.show()




##### how to transmit + other properties of SDR 

# rx_gain0 = 10
# rx_gain1 = 10



# # Example read properties
# print("RX LO %s" % (sdr.fc))

# tx_lo = fc
# tx_gain0 = -10
# tx_gain1 = -10

# # Program the Tx with some data
# fs = int(sdr.sample_rate)
# fc0 = int(2e6)
# fc1 = int(5e6)
# N = 2**16
# ts = 1 / float(fs)
# t = np.arange(0, N * ts, ts)
# i0 = np.cos(2 * np.pi * t * fc0) * 2 ** 14
# q0 = np.sin(2 * np.pi * t * fc0) * 2 ** 14
# i1 = np.cos(2 * np.pi * t * fc1) * 2 ** 14
# q1 = np.sin(2 * np.pi * t * fc1) * 2 ** 14
# iq0 = i0 + 1j * q0
# iq1 = i1 + 1j * q1
# sdr.tx([iq0, iq1])   # Send Tx data.

# from scipy import signal
# import numpy as np

# def fm_demod(rx_data, fs, fc):
    # # Channel filter
    # cutoffs = np.array([fc-100e3, fc+100e3])/(fs/2)
    # b = signal.firwin(101, cutoffs, fs=fs, pass_zero=False)
    # filtered = signal.lfilter(b, 1.0, rx_data)
    
    # # PLL/carrier recovery
    # # pll = signal.remez(64, [0, 75e3], [1], fs=fs)
    # pll = signal.remez(64, 
    #               bands=[0, 50e3, 100e3, fs/2],  # Wider transition band
    #               desired=[1, 0],
    #               fs=fs)
    # carrier = signal.lfilter(pll, 1.0, filtered)
    
    # # FM demodulation
    # phase = np.unwrap(np.angle(carrier))
    # demod = np.diff(phase)/(2*np.pi)
    
    # # Audio filter
    # audio_lpf = signal.firwin(101, 15e3, fs=fs)
    # audio = signal.lfilter(audio_lpf, 1.0, demod)
    
    # return audio, carrier
# import sk_dsp_comm 
# from commpy.equalizers import CMAEqualizer
# eq = CMAEqualizer(num_taps=64, step_size=0.01)
# cleaned_signal = eq.equalize(rxdata[0])

# def cma_clean(ref_signal, step_size=0.001, filter_len=64, iterations=100):
#     # Normalize input
#     ref_signal = ref_signal / np.sqrt(np.mean(np.abs(ref_signal)**2))
    
#     # Initialize filter weights
#     w = np.zeros(filter_len, dtype=complex)
#     w[0] = 1.0
    
#     # Create convolution matrix
#     X = np.zeros((len(ref_signal) - filter_len + 1, filter_len), dtype=complex)
#     for i in range(len(ref_signal) - filter_len + 1):
#         X[i,:] = ref_signal[i:i + filter_len]
    
#     # Track variance for debugging
#     var_history = []
    
#     # CMA iteration
#     for iter in range(iterations):
#         y = X @ w
        
#         # Modified error calculation
#         e = y * (1.0 - np.abs(y)**2)
        
#         # Update weights
#         w = w + step_size * X.conj().T @ e
#         w = w / np.linalg.norm(w)
        
#         # Track progress
#         var_history.append(np.var(np.abs(y)))
        
#     cleaned_signal = np.convolve(ref_signal, w, mode='valid')
    
#     # # Plot convergence
#     # plt.figure()
#     # plt.plot(var_history)
#     # plt.title('Modulus Variance Over Iterations')
#     # plt.show()
    
#     return cleaned_signal, w
# # Use it with your data:
# dirty_ref = rxdata[0][0:1000]
# cleaned_ref = cma_clean(dirty_ref, step_size=0.001,iterations=5)

# # Plot before/after
# fig = make_subplots(rows=1, cols=1, 
#                     subplot_titles=('Original Reference Signal', 'CMA Cleaned Signal'))
# fig.add_trace(go.Scatter(y=real(dirty_ref)), row=1, col=1)
# fig.add_trace(go.Scatter(y=real(cleaned_ref)), row=1, col=1)
# fig.update_layout(height=800, title_text="CMA Reference Signal Cleaning")
# fig.show()

# # Check modulus variance (should be smaller after cleaning)
# print(f"Original modulus variance: {np.var(np.abs(dirty_ref))}")
# print(f"Cleaned modulus variance: {np.var(np.abs(cleaned_ref))}")



# from gnuradio import gr, analog, filter, blocks
# import numpy as np

# class FMReceiver(gr.top_block):
#     def __init__(self, sample_rate=1e6, center_freq=100e6):
#         gr.top_block.__init__(self)
        
#         # Parameters
#         audio_rate = 48000
#         fm_deviation = 75e3
#         channel_width = 200e3
        
#         # Create source from your Pluto data
#         self.source = blocks.vector_source_c(rxdata[0])
        
#         # Low pass filter to isolate channel
#         self.channel_filter = filter.fir_filter_ccf(
#             1,
#             filter.firdes.low_pass(
#                 1,
#                 sample_rate,
#                 channel_width/2,
#                 channel_width/10
#             )
#         )
        
#         # PLL for carrier recovery
#         # Parameters: loop_bw, max_freq, min_freq
#         self.pll = analog.pll_carriertracking_cc(0.1, 2*np.pi*100e3, -2*np.pi*100e3)
        
#         # FM demodulator
#         self.fm_demod = analog.quadrature_demod_cf(sample_rate/(2*np.pi*fm_deviation))
        
#         # Audio filter
#         audio_decim = int(sample_rate/audio_rate)
#         self.audio_filter = filter.fir_filter_fff(
#             audio_decim,
#             filter.firdes.low_pass(
#                 1,
#                 sample_rate,
#                 15e3,
#                 1e3
#             )
#         )
        
#         # Connect blocks
#         self.connect(self.source, self.channel_filter)
#         self.connect(self.channel_filter, self.pll)
#         self.connect(self.pll, self.fm_demod)
#         self.connect(self.fm_demod, self.audio_filter)

# # Create and run flow graph
# fg = FMReceiver(sample_rate=int(fs), center_freq=int(fc))
# fg.start()

# # Get output by adding a sink:
# sink = blocks.vector_sink_f()
# fg.connect(fg.audio_filter, sink)
# fg.wait()
# demodulated = np.array(sink.data())

# # Plot results
# fig = go.Figure()
# fig.add_trace(go.Scatter(y=demodulated))
# fig.update_layout(title="Demodulated FM Signal")
# fig.show()

# # You can also look at the PLL output to check carrier recovery
# pll_sink = blocks.vector_sink_c()
# fg.connect(fg.pll, pll_sink)
# pll_out = np.array(pll_sink.data())

# # Plot PLL phase tracking
# fig = make_subplots(rows=2, cols=1, 
#                     subplot_titles=('PLL Phase', 'PLL Frequency Error'))
# fig.add_trace(go.Scatter(y=np.angle(pll_out)), row=1, col=1)
# fig.add_trace(go.Scatter(y=np.diff(np.unwrap(np.angle(pll_out)))), row=2, col=1)
# fig.update_layout(height=800, title_text="PLL Performance")
# fig.show()
