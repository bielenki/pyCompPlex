import numpy as np
from numpy.lib.stride_tricks import as_strided
import numba
from numba import jit
from osgeo import gdal
import multiprocessing as mp

#if __name__ == '__main__':
@jit
def pad_with(vector, pad_width, iaxis, kwargs):
        pad_value = kwargs.get('padder', 0)
        vector[:pad_width[0]] = pad_value
        vector[-pad_width[1]:] = pad_value
@jit
def expandStrider (input,k):
    kernel = np.arange(k*k).reshape(k, k)
    expanded_input = as_strided(
        input,
        shape=(
            input.shape[0] - kernel.shape[0] + 1,  # The feature map is a few pixels smaller than the input
            input.shape[1] - kernel.shape[1] + 1,
            kernel.shape[0],
            kernel.shape[1],
        ),
        strides=(
            input.strides[0],
            input.strides[1],
            input.strides[0],  # When we move one step in the 3rd dimension, we should move one step in the original data too
            input.strides[1],
        ),
        writeable=False,  # totally use this to avoid writing to memory in weird places
    )
    return expanded_input
@jit
def calcProb (X):
    lenVet = X.size
    prob = [(X[X==i]).size/(lenVet*1.0) for i in set(X)]
    return np.array(prob)
@jit
def convolucaoP (p):
    Hmax = np.log2(len(p)*1.0)
    ent =-p * np.log2(p)
    sEnt = sum(ent)
    C = sEnt/Hmax if (Hmax != 0.0) else 0.0
    D = sum ((p-(1/(len(p))))**2)
    SDL=(1-C)*C
    LMC=D*C
    return (np.array(sEnt), np.array(C), np.array(SDL), np.array(LMC))
