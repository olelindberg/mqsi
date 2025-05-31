import numpy as np
from mvc_integrand   import mvc_integrand
from hermite_quintic import hermite_quintic
from curve_metrics   import arc_length


def mvc_objective_function(x):
    
    ds = arc_length(x)



    gauss_xi, gauss_w = np.polynomial.legendre.leggauss(20)
    t   = (gauss_xi+1)/2
    w   = 0.5*gauss_w
    Ht   = hermite_quintic(1, t)
    Htt  = hermite_quintic(2, t)
    Httt = hermite_quintic(3, t)

    f = 0.0
    for i in range(0,len(x)-6,6):

        ii = int(i/6)

        cx = np.zeros(6)
        cx[0:3] = x[i+0:i+3]
        cx[3:6] = x[i+6:i+9]

        cy = np.zeros(6)
        cy[0:3] = x[i+3:i+6]
        cy[3:6] = x[i+9:i+12]

        cx[1] = ds[ii]*cx[1]
        cx[4] = ds[ii]*cx[4]
        cx[2] = ds[ii]**2*cx[2]
        cx[5] = ds[ii]**2*cx[5]

        cy[1] = ds[ii]*cy[1]
        cy[4] = ds[ii]*cy[4]
        cx[2] = ds[ii]**2*cy[2]
        cx[5] = ds[ii]**2*cy[5]

        cx_t   = Ht@cx
        cy_t   = Ht@cy
        cx_tt  = Htt@cx
        cy_tt  = Htt@cy
        cx_ttt = Httt@cx
        cy_ttt = Httt@cy

        integrand = mvc_integrand(cx_t, cy_t, cx_tt, cy_tt, cx_ttt, cy_ttt)    

        f = f + np.dot(integrand,w)

    return f
