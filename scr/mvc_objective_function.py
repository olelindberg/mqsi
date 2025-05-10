import numpy as np
from mvc_integrand   import mvc_integrand
from mvc_vertex_to_curve import mvc_vertex_to_curve
from hermite_quintic import hermite_quintic

def mvc_objective_function(x,s):

    gauss_xi, gauss_w = np.polynomial.legendre.leggauss(20)


    f = 0.0
    for i in range(0,len(x)-6,6):

        x = np.reshape(x,(int(len(x)/6),6))
        ii = int(i/6)
        cx,cy = mvc_vertex_to_curve(x,ii)
        x = x.flatten()

        s_xi = 0.5*(s[ii+1] - s[ii])
        Ht   = hermite_quintic(1, gauss_xi,s_xi)
        Htt  = hermite_quintic(2, gauss_xi,s_xi)
        Httt = hermite_quintic(3, gauss_xi,s_xi)


        cx_t   = Ht@cx
        cy_t   = Ht@cy
        cx_tt  = Htt@cx
        cy_tt  = Htt@cy
        cx_ttt = Httt@cx
        cy_ttt = Httt@cy

        integrand = mvc_integrand(cx_t, cy_t, cx_tt, cy_tt, cx_ttt, cy_ttt)    

        f = f + np.sum(integrand*gauss_w)

    return f
