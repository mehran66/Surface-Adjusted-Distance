import itertools
import numpy as np
import matplotlib.pyplot as plt

# This module can be used to fit a polynomial surface to a set of (x,y,z) points
#Reference: http://stackoverflow.com/questions/7997152/python-3d-polynomial-surface-fit-order-dependent
# polyfit2d funcion is used first to find the coefficient of a polynomila of order m fitted to known coordinates x, y and z
# polyval2d is used to interpolate the coordinate of an unknown point

def polyfit2d(x, y, z, order=3):
    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, _, _, _ = np.linalg.lstsq(G, z)
    return m

def polyval2d(x, y, m):
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z




