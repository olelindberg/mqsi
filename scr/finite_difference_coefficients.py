import numpy as np
import math

def finite_difference_coefficients(x,x0,derivative):

    n = len(x)
    M = np.zeros((n,n))

    for i in range(n):
        for j in range(n):
        
            M[i,j] = 1.0/math.factorial(j)*(x[i]-x0)**j

    M_inv = np.linalg.inv(M)
    
    return M_inv[derivative,:]
