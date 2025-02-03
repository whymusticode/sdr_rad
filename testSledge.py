

 





import numpy as np
import os
import time
import importlib

mydir = '/mnt/windows/BACKUP/CODE/SDR'
# mydir='/home/bill/Desktop'
os.chdir(mydir)

import sledge
# import ECA 

p = dict()
p['fs'] = 2e5
tau = 3
p['fc'] = 1e8
p['vMax'] = 400
p['rMax'] = 300e3
N = int(round(2e5*tau))
Ntaps=30
tx = np.random.randn(N) + 1j*np.random.randn(N) 

delay = 100
rx = np.concatenate([np.zeros(delay),tx[:-delay]],axis=0)

delay2 = 50
rx2 = np.concatenate([np.zeros(delay2),tx[:-delay2]],axis=0)

rx+=rx2*1e-2



# taps = np.random.randn(Ntaps)
# tmp =  np.convolve(tx, taps, mode='full')*1e3
# rx += tx*1e5

# rx = rx[40:-40]
# importlib.reload(ECA)
# start = time.time()
# [rxClean, tapEst, cancelDB] = ECA.ECA(tx, rx[40:-40],tapEst0=tapEst)
# print(time.time()-start)
# importlib.reload(sledge)
# rx = np.concatenate((np.zeros(40),rxClean,np.zeros(40)),axis=0)

start = time.time()
(range_vals, velocity, image) = sledge.sledge(tx, rx, p)
end = time.time()
print("Elapsed time: ", end - start)

image.shape


z = 10*np.log10(image.real**2 + image.imag**2)

import plotly.graph_objects as go
fig = go.Figure(data=go.Heatmap(z=z,x=range_vals/1e3,y=velocity,))
fig.update_layout(template='plotly_dark',hovermode=False)
fig.show()

# device = np.device("cuda" if np.cuda.is_available() else "cpu")
# device = "cpu"
# device = 'cuda'
# rx2 = np.cat([np.zeros(delay2),tx[:-delay2]],dim=0)
# np.cuda.synchronize()


# N = 100
# delay = 50

# tx = np.randn(N) + 1j*np.randn(N) 
# rx = np.cat([np.zeros(delay),tx[:-delay]],dim=0)



# import numpy as np

# def create_circshift_matrix(data, num_rows):
#     data = np.array(data)
#     matrix = np.vstack([np.roll(data, -i) for i in range(num_rows)])
#     return matrix

# result = create_circshift_matrix(tx, N)
# result = np.tensor(result[2:,:])
# new = np.conj(result.T) @ result
# txMod = new**-1 @ tx

# y = np.fft.ifft(np.conj(np.fft.fft(tx))*np.fft.fft(rx))
# y2 = np.fft.ifft(np.conj(np.fft.fft(txMod))*np.fft.fft(rx))

# ydb = 10*np.log10(y.real**2 + y.imag**2)
# ydb2 = 10*np.log10(y2.real**2 + y2.imag**2)

# # Create a figure with two scatter plots
# fig = go.Figure(
#     data=[
#         go.Scatter(
#             x=np.arange(0, N).numpy(),
#             y=ydb.cpu(),
#             name='Scatter 1'
#         ),
#         go.Scatter(
#             x=np.arange(0, N).numpy(),
#             y=ydb2.cpu(),
#             name='Scatter 2'
#         )
#     ]
# )
# fig.update_layout(template='plotly_dark',hovermode=False)

# fig.show()






# import matplotlib.pyplot as plt

# plt.imshow(imdb.cpu(), aspect='auto', 
#            extent=[range_vals.min()/1e3, range_vals.max()/1e3, velocity.min(), velocity.max()])
# plt.xlabel('range (km)')
# plt.ylabel('velocity (m/s)')
# plt.colorbar()
# plt.show()


# imdb_list_of_vectors = imdb.cpu().numpy().tolist()

# imdb_list_of_vectors = imdb.cpu().numpy().T.tolist()


# import plotly.graph_objects as go

# fig = px.imshow(
#     imdb_list_of_vectors,
#     labels=dict(x="range (km)", y="velocity (m/s)"),
#     x = range_vals/1e3,
#     y = velocity,
#     template='plotly_dark',
# )
# fig.update_layout(
#     xaxis=dict(scaleanchor=None),
#     yaxis=dict( scaleanchor=None ),
# )
# fig.show()