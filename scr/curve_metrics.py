import numpy as np
from mvc_vertex_to_curve import mvc_vertex_to_curve
from hermite_quintic_derivatives import hermite_quintic_derivatives

def tangent_vector_length(x_s,y_s):
    return np.sqrt(x_s**2 + y_s**2)

def curvature(x_ss,y_ss):
    return np.sqrt(x_ss**2 + y_ss**2)

def arc_length(x,debug=False):
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

        ds=1
        for k in range(100):
            
            x = np.reshape(x,(num_points,6))
            ii = int(i/6)
            cx,cy = mvc_vertex_to_curve(x,ii)
            x = x.flatten()

            cx[1] = ds*cx[1]
            cx[4] = ds*cx[4]
            cx[2] = ds**2*cx[2]
            cx[5] = ds**2*cx[5]

            cy[1] = ds*cy[1]
            cy[4] = ds*cy[4]
            cy[2] = ds**2*cy[2]
            cy[5] = ds**2*cy[5]

            xx, yy, x_t, y_t, x_tt, y_tt,x_ttt,y_ttt = hermite_quintic_derivatives(1,t,cx,cy)

            ds = np.dot(np.sqrt(x_t**2 + y_t**2),w)
            if (debug):
                print("ds = ", ds)
        l.append(ds)


    return l