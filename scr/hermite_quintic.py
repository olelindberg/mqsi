import numpy as np

def hermite_quintic(nderivative, x):
    """
    Hermite quintic polynomial basis functions and their derivatives.
    The polynomials are defined in the range [-1, 1]:

    p = C0*H0 + C0t*H1 + C0tt*H2 + C1*H3 + C1t*H4 + C1tt*H5
    
    Parameters
    ----------
    nderivative : int
        The order of the derivative to compute.
    x : array_like
        The input values at which to evaluate the basis functions.
    Returns
    -------
    H : array_like
        The evaluated basis functions or their derivatives.
    """
    x = np.asarray(x)
    H = np.zeros((6, len(x)))

    # Hermite quintic polynomial coefficients
    if nderivative == 0:
        H[0, :] = 1/2 - (15*x)/16 + (5*x**3)/8 - (3*x**5)/16
        H[1, :] = 5/16 - (7*x)/16 - (3*x**2)/8 + (5*x**3)/8 + x**4/16 - (3*x**5)/16
        H[2, :] = 1/16 - x/16 - x**2/8 + x**3/8 + x**4/16 - x**5/16
        H[3, :] = 1/2 + (15*x)/16 - (5*x**3)/8 + (3*x**5)/16
        H[4, :] = -5/16 - (7*x)/16 + (3*x**2)/8 + (5*x**3)/8 - x**4/16 - (3*x**5)/16
        H[5, :] = 1/16 + x/16 - x**2/8 - x**3/8 + x**4/16 + x**5/16

    elif nderivative == 1:
        H[0, :] = -15/16 + (15*x**2)/8 - (15*x**4)/16
        H[1, :] = -7/16 - (3*x)/4 + (15*x**2)/8 + x**3/4 - (15*x**4)/16
        H[2, :] = -1/16 - x/4 + (3*x**2)/8 + x**3/4 - (5*x**4)/16
        H[3, :] = 15/16 - (15*x**2)/8 + (15*x**4)/16
        H[4, :] = -7/16 + (3*x)/4 + (15*x**2)/8 - x**3/4 - (15*x**4)/16
        H[5, :] = 1/16 - x/4 - (3*x**2)/8 + x**3/4 + (5*x**4)/16

    elif nderivative == 2:
        H[0, :] = (15*x)/4 - (15*x**3)/4
        H[1, :] = -3/4 + (15*x)/4 + (3*x**2)/4 - (15*x**3)/4
        H[2, :] = -1/4 + (3*x)/4 + (3*x**2)/4 - (5*x**3)/4
        H[3, :] = -15*x/4 + (15*x**3)/4
        H[4, :] = 3/4 + (15*x)/4 - (3*x**2)/4 - (15*x**3)/4
        H[5, :] = -1/4 - (3*x)/4 + (3*x**2)/4 + (5*x**3)/4

    elif nderivative == 3:
        H[0, :] = 15/4 - (45*x**2)/4
        H[1, :] = 15/4 + (6*x)/4 - (45*x**2)/4
        H[2, :] = 3/4 + (6*x)/4 - (15*x**2)/4
        H[3, :] = -15/4 + (45*x**2)/4
        H[4, :] = 15/4 - (6*x)/4 - (45*x**2)/4
        H[5, :] = -3/4 + (6*x)/4 + (15*x**2)/4

    # Transpose to match MATLAB's H = H'
    return H.T