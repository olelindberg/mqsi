import numpy as np

def hermite_quintic(nderivative, t,dsdt=1):
    """
    Hermite quintic polynomial basis functions and their derivatives.
    The polynomials are defined in the range [-1, 1]:

    p = C0*H0 + C0t*H1 + C0tt*H2 + C1*H3 + C1t*H4 + C1tt*H5
    
    Parameters
    ----------
    nderivative : int
        The order of the derivative to compute.
    t : array_like
        The input values at which to evaluate the basis functions.
    Returns
    -------
    H : array_like
        The evaluated basis functions or their derivatives.
    """
    t = np.asarray(t)
    H = np.zeros((6, len(t)))

    # Hermite quintic polynomial coefficients
    if nderivative == 0:
        H[0,:] = 1 -                 10*t**3 +  15*t**4 -   6*t**5     
        H[1,:] =     t -              6*t**3 +   8*t**4 -   3*t**5       
        H[2,:] =         0.5*t**2 - 1.5*t**3 + 1.5*t**4 - 0.5*t**5 
        H[3,:] =                     10*t**3 -  15*t**4 +   6*t**5       
        H[4,:] =                  -   4*t**3 +   7*t**4 -   3*t**5        
        H[5,:] =                    0.5*t**3 -     t**4 + 0.5*t**5
    elif nderivative == 1:

        H[0,:] = -30*t**4 + 60*t**3 - 30*t**2
        H[1,:] = -15*t**4 + 32*t**3 - 18*t**2 + 1
        H[2,:] = -2.5*t**4 + 6.0*t**3 - 4.5*t**2 + 1.0*t
        H[3,:] = 30*t**4 - 60*t**3 + 30*t**2
        H[4,:] = -15*t**4 + 28*t**3 - 12*t**2
        H[5,:] = 2.5*t**4 - 4*t**3 + 1.5*t**2

    elif nderivative == 2:

        H[0,:] = -120*t**3 + 180*t**2 - 60*t
        H[1,:] = -60*t**3 + 96*t**2 - 36*t
        H[2,:] = -10.0*t**3 + 18.0*t**2 - 9.0*t + 1.0
        H[3,:] = 120*t**3 - 180*t**2 + 60*t
        H[4,:] = -60*t**3 + 84*t**2 - 24*t
        H[5,:] = 10.0*t**3 - 12*t**2 + 3.0*t

    elif nderivative == 3:

        H[0,:] = -360*t**2 + 360*t - 60
        H[1,:] = -180*t**2 + 192*t - 36
        H[2,:] = -30.0*t**2 + 36.0*t - 9.0
        H[3,:] = 360*t**2 - 360*t + 60
        H[4,:] = -180*t**2 + 168*t - 24
        H[5,:] = 30.0*t**2 - 24*t + 3.0

    H[0, :] = H[0, :]
    H[1, :] = H[1, :]*dsdt
    H[2, :] = H[2, :]*dsdt**2
    H[3, :] = H[3, :]
    H[4, :] = H[4, :]*dsdt
    H[5, :] = H[5, :]*dsdt**2

    # Transpose to match MATLAB's H = H'
    return H.T