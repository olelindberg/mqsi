import numpy as np
from mvc_integrand   import mvc_integrand
from hermite_quintic import hermite_quintic

def mvc_objective_function(x):

    gauss_xi, gauss_w = np.polynomial.legendre.leggauss(20)

    Ht   = hermite_quintic(1, gauss_xi)
    Htt  = hermite_quintic(2, gauss_xi)
    Httt = hermite_quintic(3, gauss_xi)

    f = 0.0
    for i in range(0,len(x)-6,6):

        cx = np.zeros(6)
        cx[0:3] = x[i+0:i+3]
        cx[3:6] = x[i+6:i+9]

        cy = np.zeros(6)
        cy[0:3] = x[i+3:i+6]
        cy[3:6] = x[i+9:i+12]

        cx_t   = Ht@cx
        cy_t   = Ht@cy
        cx_tt  = Htt@cx
        cy_tt  = Htt@cy
        cx_ttt = Httt@cx
        cy_ttt = Httt@cy

        integrand = mvc_integrand(cx_t, cy_t, cx_tt, cy_tt, cx_ttt, cy_ttt)    

        f = f + np.sum(integrand*gauss_w)

    return f
