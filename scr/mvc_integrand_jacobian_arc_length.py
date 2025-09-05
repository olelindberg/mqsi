import numpy as np
from mvc_integrand_gradient import mvc_integrand_gradient
from hermite_quintic        import hermite_quintic

def mvc_integrand_jacobian_arc_length(x,ds_all):

    jac  = 0*ds_all

    gauss_xi, gauss_w = np.polynomial.legendre.leggauss(20)

    t    = (gauss_xi+1)/2
    Ht   = hermite_quintic(1, t)
    Htt  = hermite_quintic(2, t)
    Httt = hermite_quintic(3, t)

    for i in range(0,len(x)-6,6):

        ii = int(i/6)

        ds = ds_all[ii]

        cx = np.zeros(6)
        cx[0:3] = x[i+0:i+3]
        cx[3:6] = x[i+6:i+9]

        cy = np.zeros(6)
        cy[0:3] = x[i+3:i+6]
        cy[3:6] = x[i+9:i+12]

        cx[1] = ds*cx[1]
        cx[4] = ds*cx[4]
        cx[2] = ds**2*cx[2]
        cx[5] = ds**2*cx[5]
        
        cy[1] = ds*cy[1]
        cy[4] = ds*cy[4]
        cy[2] = ds**2*cy[2]
        cy[5] = ds**2*cy[5]

        cx_t   = Ht@cx
        cy_t   = Ht@cy
        cx_tt  = Htt@cx
        cy_tt  = Htt@cy
        cx_ttt = Httt@cx
        cy_ttt = Httt@cy

        # Derivatives with respect to arc length

        cx = np.zeros(6)
        cx[0:3] = x[i+0:i+3]
        cx[3:6] = x[i+6:i+9]

        cy = np.zeros(6)
        cy[0:3] = x[i+3:i+6]
        cy[3:6] = x[i+9:i+12]

        cx[0] = 0
        cx[1] = cx[1]
        cx[2] = 2*ds*cx[2]
        cx[3] = 0
        cx[4] = cx[4]
        cx[5] = 2*ds*cx[5]

        cy[0] = 0
        cy[1] = cy[1]
        cy[2] = 2*ds*cy[2]
        cy[3] = 0
        cy[4] = cy[4]
        cy[5] = 2*ds*cy[5]

        cx_at   = Ht@cx
        cy_at   = Ht@cy
        cx_att  = Htt@cx
        cy_att  = Htt@cy
        cx_attt = Httt@cx
        cy_attt = Httt@cy

        jac[ii] = np.sum(mvc_integrand_gradient(cx_t,cy_t,cx_tt,cy_tt,cx_ttt,cy_ttt,cx_at,cy_at,cx_att,cy_att,cx_attt,cy_attt)*gauss_w/2)
        
    return jac
