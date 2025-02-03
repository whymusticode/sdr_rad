

import os
os.chdir('/mnt/windows/BACKUP/CODE/SDR')
# import resampler
import time
import numpy as np

# from resample import resampstr







import numpy as np
import math

def resampstr(p, m=-1):
    n = len(p)
    if m == -1:
        m = n

    # Normalize p to sum to m
    pn = np.array(p) / np.sum(p) * m

    # Initialize output and random values
    s = np.zeros(m, dtype=np.int32)
    r = np.random.rand(n)

    # Resampling logic
    c = 0.0
    k = 0
    for i in range(n):
        c += pn[i]
        if c >= 1.0:
            a = math.floor(c)
            c -= a
            for j in range(int(a)):
                s[k] = i + 1  # MATLAB uses 1-based indexing
                k += 1
        if k < m and c >= r[k]:
            c -= 1.0
            s[k] = i + 1
            k += 1

    # Handle last element
    if s[m - 1] == 0:
        mode_val = np.argmax(np.bincount(s))
        s[m - 1] = mode_val

    return s




m = 10000000
p = np.random.rand(m)

# p = [0.2, 0.3, 0.5]


# Call the C++ function from Python
start = time.time()
s =resampstr(p, m)
Tused = time.time()-start
print(str(Tused)+' secs 2 run')

m/Tused/1e6
# print(s)