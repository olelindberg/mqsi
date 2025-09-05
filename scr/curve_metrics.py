import numpy as np
from mvc_vertex_to_curve import mvc_vertex_to_curve
from hermite_quintic_derivatives import hermite_quintic_derivatives
from mvc_integrand_jacobian_arc_length import mvc_integrand_jacobian_arc_length

def tangent_vector_length(x_s,y_s):
    return np.sqrt(x_s**2 + y_s**2)

def curvature(x_ss,y_ss):
    return np.sqrt(x_ss**2 + y_ss**2)

def arc_length(x,ds,debug=False,tol=1e-15):
    """
    xi is the gauss quadrature points
    t  is the hermite quintic parameter 
    """

    if debug:
        print("Computing arc length ...")
    
    gausss_xi, gauss_w = np.polynomial.legendre.leggauss(24)
    t = 0.5*(gausss_xi+1)
    w = 0.5*gauss_w

    l = []

    num_points = int(len(x)/6)



    ds_old = 0*ds
    grad   = 0*ds
    for k in range(100):


        grad_old = grad

        grad     = mvc_integrand_jacobian_arc_length(x,ds)
        
        #----------------------------------------------------#
        # Compute iteration step size:
        #----------------------------------------------------#
        dds = ds - ds_old
        dgrad = grad - grad_old
#        print(dgrad)
        gamma = 1e-10
        if dgrad.dot(dgrad)>0:
            gamma = np.abs(dds.dot(dgrad))/dgrad.dot(dgrad)

        #----------------------------------------------------#
        # Update the solution:
        #----------------------------------------------------#
        ds_old = ds
        ds     = ds + gamma*grad
#        print(f"  Iter {k:<3} gamma: {gamma:<10.4e}")
#        print(grad)
        print(ds)

    return ds