
# kernprof -l -v sledge.py
# python3 sledge.py

# import torch
import numpy as np

# from numba import jit

from scipy.constants import c

import math

def nearest_power_of_2(x):
    if x <= 0:
        raise ValueError("Input must be a positive integer")
    return 1 << (x - 1).bit_length()

# @jit(nopython=True, parallel=True)
def sledge(tx, rx, p):
    """
    Performs 2-D Cross Correlation Function in range and Doppler for
    arbitrary transmit waveform. Caution, this has Doppler ambiguities.

    Inputs:
        tx - IQ signal transmitted (numpy array)
        rx - IQ signal received (numpy array), length(rx) == length(tx)
        p - parameters including:
            p['fs']: sample frequency (Hz)
            p['fc']: center frequency (Hz)
            p['rMax']: maximum range to create image (m)
            p['vMax']: max velocity +- to create image (m/s)

    Outputs:
        range_vals - (m)
        velocity - (m/s)
        image - complex range-Doppler image, [velocity x range]
    """
    # Constants
    n0 = len(tx)  # length of inputs
    dpp = c / n0 * p['fs'] / p['fc']  # Doppler per pixel approximation
    nD = round(2 * p['vMax'] / dpp)  # number of Doppler bins
    # nD = nearest_power_of_2(nD)
    nB = n0 // nD  # batch length
    n = nD * nB  # new integration length
    dpp = c / n * p['fs'] / p['fc']  # Doppler per pixel exact
    rpp = c / p['fs']  # range per pixel
    range_vals = np.arange(rpp, p['rMax'] + rpp, rpp)  # range bins
    nR = len(range_vals)  # number of range bins
    velocity = np.fft.fftshift(np.fft.fftfreq(nD, d=1/nD)) * -dpp
    nBR = nB + nR  # length of convolution
    

    part1 = rx[:n].reshape(nD, nB)
    part2 = np.concatenate([part1[1:, :nR], np.zeros((1, nR), dtype=rx.dtype)], axis=0)
    rxBatch = np.concatenate([part1, part2], axis=1)

    TX = np.conj(np.fft.fft(tx[:n].reshape(nD, nB), n=nBR, axis=1))
    RX = np.fft.fft(rxBatch, axis=1)
    corFun1d = np.fft.ifft(TX * RX, axis=1)
    image = np.fft.fftshift(np.fft.fft(corFun1d[:, :nR], axis=0), axes=0)

    return range_vals, velocity, image

# # @profile
# def sledge(tx, rx, p):
#     """
#     Performs 2-D Cross Correlation Function in range and Doppler for
#     arbitrary transmit waveform. Caution, this has Doppler ambiguities.

#     Inputs:
#         tx - IQ signal transmitted (Tensor)
#         rx - IQ signal received (Tensor), length(rx) == length(tx)
#         p - parameters including:
#             p['fs']: sample frequency (Hz)
#             p['fc']: center frequency (Hz)
#             p['rMax']: maximum range to create image (m)
#             p['vMax']: max velocity +- to create image (m/s)

#     Outputs:
#         range - (m)
#         velocity - (m/s)
#         image - complex range-Doppler image, [velocity x range]
#     """
#     # Constants
#     n0 = len(tx)  # length of inputs
#     dpp = c / n0 * p['fs'] / p['fc']  # Doppler per pixel approximation
#     nD = round(2 * p['vMax'] / dpp)  # number of Doppler bins
#     nB = n0 // nD  # batch length
#     n = nD * nB  # new integration length
#     dpp = c / n * p['fs'] / p['fc']  # Doppler per pixel exact
#     rpp = c / p['fs']  # range per pixel
#     range_vals = torch.arange(rpp, p['rMax'] + rpp, rpp)  # range bins
#     nR = len(range_vals)  # number of range bins
#     velocity = torch.fft.fftshift(torch.fft.fftfreq(nD,d=1/nD))* -dpp
#     nBR = nB + nR  # length of convolution
    
#     part1 = rx[:n].view(nD,nB)
#     part2 = torch.cat([part1[1:,:nR],torch.zeros(1,nR,device = rx.device,dtype=rx.dtype)], dim=0)
#     rxBatch = torch.cat([part1, part2], dim=1)

#     TX = torch.conj(torch.fft.fft(tx[:n].view(nD, nB) , n=nBR, dim=1) )
#     RX = torch.fft.fft(rxBatch, dim=1)
#     corFun1d = torch.fft.ifft(TX*RX, dim=1)
#     image = torch.fft.fftshift(torch.fft.fft(corFun1d[:,:nR], dim=0), dim=0)

#     return range_vals, velocity, image


# nD = 10
# torch.arange(-nD // 2, nD // 2)
# frequencies = torch.fft.ifftshift(torch.fft.fftfreq(nD,d=1/nD))* -dpp



# fig.update_layout(
#     # xaxis_title='range (km)',     # Uncomment and modify as needed
#     # yaxis_title='velocity (mps)', # Uncomment and modify as needed
#     # coloraxis_colorbar=dict(title='Intensity')  # Example colorbar label, can be modified
# )





# image = torch.fft.ifft(   torch.conj(torch.fft.fft(tx))   * torch.fft.fft(rx)   )
# z = image.real**2 + image.imag**2 
# imdb = 10*torch.log10(z)
# fig = px.scatter(x=torch.arange(0,N), y=imdb.cpu(), title="Quick Scatter Plot")
# fig.show()
# image = corFun1d[:,:nR]
# image = torch.fft.fft(corFun1d[:,:nR], dim=0)
# image = torch.fft.fftshift(torch.fft.fft(corFun1d[:,:nR], dim=0), dim=0)


