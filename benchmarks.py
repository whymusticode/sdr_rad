
import multiprocessing

multiprocessing.cpu_count()


import scipy
import numpy as np
from scipy.fft import fft
import time
# from scipy.stats import norm
# x = norm.rvs(size=sz,) + 1j*norm.rvs(size=sz,)
import torch



ct = 0
sz = (256,1028)
x = np.random.random(sz) + 1j*np.random.random(sz,)
x = torch.tensor(x,dtype=torch.complex64)
start = time.time()
while time.time()-start < 10:
    ct += 1
    
    # x = norm.rvs(size=sz,) + 1j*norm.rvs(size=sz,)

    # y = fft(x,axis=0) #,overwrite_x=True
    # y = np.fft.fft(x,axis=0)
    y = torch.fft.fft(x,dim=0)
print(ct)







# # Example usage and benchmark
# if __name__ == '__main__':
#     import time
    
#     # Test data
#     M, N = 1000, 2000  # Adjust sizes based on your needs
#     test_data = np.random.random((M, N)) + 1j * np.random.random((M, N))
    
#     # For 4-core ARM, you might want to explicitly specify CPU IDs
#     # Example: Using cores 0,1,2,3
#     cpu_ids = [0, 1, 2, 3]
    
#     print("Starting benchmark...")
    
#     # Serial FFT
#     start = time.time()
#     serial_result = np.fft.fft(test_data, axis=0)
#     serial_time = time.time() - start
#     print(f"Serial FFT time: {serial_time:.3f} seconds")
    
#     # Parallel FFT with affinity
#     start = time.time()
#     parallel_result = parallel_fft_with_affinity(test_data, n_threads=4, cpu_ids=cpu_ids)
#     parallel_time = time.time() - start
#     print(f"Parallel FFT time with affinity: {parallel_time:.3f} seconds")
    
#     # Verify results
#     max_diff = np.max(np.abs(serial_result - parallel_result))
#     print(f"Maximum difference between results: {max_diff}")
#     print(f"Speedup: {serial_time/parallel_time:.2f}x")

#     # Show current system CPU info
#     if psutil is not None:
#         print("\nCPU Information:")
#         print(f"Physical cores: {psutil.cpu_count(logical=False)}")
#         print(f"Total cores: {psutil.cpu_count(logical=True)}")
#         print(f"Current CPU frequencies:")
#         freq = psutil.cpu_freq(percpu=True)
#         if freq:
#             for i, f in enumerate(freq):
#                 print(f"CPU {i}: {f.current:.0f} MHz")


