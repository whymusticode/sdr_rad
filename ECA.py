import numpy as np
from scipy.sparse.linalg import LinearOperator, lsqr

def ECA(tx, rx, tapEst0=None):
    """
    Solves for the TAPS that minimize ||rx - conv(tx,tapEst,'valid')||
    and returns rxClean = rx - conv(tx,tapEst,'valid')

    Parameters:
        tx (numpy.ndarray): Transmit signal
        rx (numpy.ndarray): Receive signal
        tapEst0 (numpy.ndarray, optional): Initial estimate for taps
    
    Returns:
        rxClean (numpy.ndarray): Cleaned receive signal
        tapEst (numpy.ndarray): Estimated taps
        cancelDB (float, optional): Cancellation in dB
    """
    # Convergence parameters
    maxIter = 50
    tolerance = 1e-3

    N = len(tx)
    Ntaps = N - len(rx) + 1
    Nfft1 = N + Ntaps - 1
    Nfft2 = N + Nfft1 - 1

    pre1 = np.fft.fft(tx, Nfft1)
    pre2 = np.fft.fft(np.flipud(np.conj(tx)), Nfft2)

    # Define the function for lsqr solver
    def ax(x, opt, pre1, pre2, N, Nfft1, Nfft2, Ntaps):
        if opt == 'notransp':
            y = np.fft.ifft(pre1 * np.fft.fft(x, Nfft1))
            y = y[Ntaps-1:len(y)-Ntaps+1]
        else:
            y = np.fft.ifft(pre2 * np.fft.fft(x, Nfft2))
            y = y[N-Ntaps:N]
        return y

    class ECALinearOperator(LinearOperator):
        def __init__(self, shape, dtype, pre1, pre2, N, Nfft1, Nfft2, Ntaps):
            super().__init__(dtype=dtype, shape=shape)
            self.pre1 = pre1
            self.pre2 = pre2
            self.N = N
            self.Nfft1 = Nfft1
            self.Nfft2 = Nfft2
            self.Ntaps = Ntaps

        def _matvec(self, x):
            return ax(x, 'notransp', self.pre1, self.pre2, self.N, self.Nfft1, self.Nfft2, self.Ntaps)

        def _rmatvec(self, x):
            return ax(x, 'transp', self.pre1, self.pre2, self.N, self.Nfft1, self.Nfft2, self.Ntaps)

    # Create a linear operator using the matrix-vector product functions
    A = ECALinearOperator((len(rx), Ntaps), dtype=np.complex128, pre1=pre1, pre2=pre2, N=N, Nfft1=Nfft1, Nfft2=Nfft2, Ntaps=Ntaps)

    # Solve for the taps
    if tapEst0 is None:
        result = lsqr(A, rx, atol=tolerance, btol=tolerance, iter_lim=maxIter)
    else:
        result = lsqr(A, rx, atol=tolerance, btol=tolerance, iter_lim=maxIter, x0=tapEst0)
    
    tapEst = result[0]
    
    # Compute clutter and cleaned signal
    rxClutter = A @ tapEst
    rxClean = rx - rxClutter

    #     # Optionally return cancellation in dB
    #     if nargout() == 3:
    cancelDB = 10 * np.log10(np.sum(np.abs(rxClean) ** 2) / np.sum(np.abs(rx) ** 2))
    print(f'Cancellation of {round(cancelDB)} dB')
    return rxClean, tapEst, cancelDB
    




#     return rxClean, tapEst

# # Helper function for nargout in Python
# def nargout():
#     import inspect
#     return len(inspect.signature(inspect.stack()[1][0].f_globals[inspect.stack()[1][0].f_code.co_name]).parameters)




# import numpy as np
# from scipy.sparse.linalg import lsqr
# from scipy.sparse.linalg import LinearOperator





# def mv(v):
#     return np.array([2*v[0], 3*v[1]])


# A = LinearOperator((2,2), matvec=mv)

# A

# A.matvec(np.ones(2))
# array([ 2.,  3.])

# A * np.ones(2)
# array([ 2.,  3.])


# def ECA(tx, rx, tapEst0=None):
#     """
#     Solves for the TAPS that minimize ||rx - conv(tx,tapEst,'valid')||
#     and returns rxClean = rx - conv(tx,tapEst,'valid')

#     Parameters:
#         tx (numpy.ndarray): Transmit signal
#         rx (numpy.ndarray): Receive signal
#         tapEst0 (numpy.ndarray, optional): Initial estimate for taps
    
#     Returns:
#         rxClean (numpy.ndarray): Cleaned receive signal
#         tapEst (numpy.ndarray): Estimated taps
#         cancelDB (float, optional): Cancellation in dB
#     """
#     # Convergence parameters
#     maxIter = 50
#     tolerance = 1e-3

#     N = len(tx)
#     Ntaps = N - len(rx) + 1
#     Nfft1 = N + Ntaps - 1
#     Nfft2 = N + Nfft1 - 1

#     pre1 = np.fft.fft(tx, Nfft1)
#     pre2 = np.fft.fft(np.flipud(np.conj(tx)), Nfft2)

#     # Define the function for lsqr solver
#     def ax(x, opt, pre1, pre2, N, Nfft1, Nfft2, Ntaps):
#         if opt == 'notransp':
#             y = np.fft.ifft(pre1 * np.fft.fft(x, Nfft1))
#             y = y[Ntaps-1:len(y)-Ntaps+1]
#         else:
#             y = np.fft.ifft(pre2 * np.fft.fft(x, Nfft2))
#             y = y[N-Ntaps:N]
#         return y

#     fun = lambda x, opt: ax(x, opt, pre1, pre2, N, Nfft1, Nfft2, Ntaps)

#     # Solve for the taps
#     if tapEst0 is None:
#         result = lsqr(lambda x: fun(x, 'notransp'), rx, atol=tolerance, btol=tolerance, iter_lim=maxIter)
#     else:
#         result = lsqr(lambda x: fun(x, 'notransp'), rx, atol=tolerance, btol=tolerance, iter_lim=maxIter, x0=tapEst0)
    
#     tapEst = result[0]
    
#     # Compute clutter and cleaned signal
#     rxClutter = fun(tapEst, 'notransp')
#     rxClean = rx - rxClutter

#     # Optionally return cancellation in dB
#     if nargout == 3:
#         cancelDB = 10 * np.log10(np.sum(np.abs(rxClean) ** 2) / np.sum(np.abs(rx) ** 2))
#         print(f'Cancellation of {round(cancelDB)} dB')
#         return rxClean, tapEst, cancelDB
    
#     return rxClean, tapEst

# # Helper function for nargout in Python
# def nargout():
#     import inspect
#     return len(inspect.signature(inspect.stack()[1][0].f_globals[inspect.stack()[1][0].f_code.co_name]).parameters)
