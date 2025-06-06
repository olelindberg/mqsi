import numpy as np
from mvc_vertex_to_curve import mvc_vertex_to_curve
from hermite_quintic_derivatives import hermite_quintic_derivatives

def tangent_vector_length(x_s,y_s):
    return np.sqrt(x_s**2 + y_s**2)

def curvature(x_ss,y_ss):
    return np.sqrt(x_ss**2 + y_ss**2)

def arc_length(x,debug=False,tol=1e-15):
    """
    xi is the gauss quadrature points
    t  is the hermite quintic parameter 
    """
    gausss_xi, gauss_w = np.polynomial.legendre.leggauss(24)
    t = 0.5*(gausss_xi+1)
    w = 0.5*gauss_w

    l = []

    num_points = int(len(x)/6)
    for i in range(0,len(x)-6,6):

        x = np.reshape(x,(num_points,6))
        ii = int(i/6)
        cx,cy = mvc_vertex_to_curve(x,ii)
        x = x.flatten()

        dx = cx[3] - cx[0]
        dy = cy[3] - cy[0]
        ds = np.sqrt(dx**2 + dy**2)
        ds_init = ds
        for k in range(100):
            
            ccx = np.zeros(6)
            ccx[0] =       cx[0]
            ccx[1] =    ds*cx[1]
            ccx[2] = ds**2*cx[2]
            ccx[3] =       cx[3]
            ccx[4] =    ds*cx[4]
            ccx[5] = ds**2*cx[5]

            ccy = np.zeros(6)
            ccy[0] =       cy[0]
            ccy[1] =    ds*cy[1]
            ccy[2] = ds**2*cy[2]
            ccy[3] =       cy[3]
            ccy[4] =    ds*cy[4]
            ccy[5] = ds**2*cy[5]

            xx, yy, x_t, y_t, x_tt, y_tt,x_ttt,y_ttt = hermite_quintic_derivatives(1,t,ccx,ccy)

            ds_old = ds
            ds = np.dot(np.sqrt(x_t**2 + y_t**2),w)
            dds_rel = np.abs(ds - ds_old)/ds_init

            if (dds_rel<tol):
                break

        l.append(ds)


    return l