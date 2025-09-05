import numpy as np
from mvc_vertex_to_curve import mvc_vertex_to_curve
from hermite_quintic_derivatives import hermite_quintic_derivatives
from mvc_integrand_jacobian_arc_length import mvc_integrand_jacobian_arc_length

def tangent_vector_length(x_s,y_s):
    return np.sqrt(x_s**2 + y_s**2)

def curvature(x_ss,y_ss):
    return np.sqrt(x_ss**2 + y_ss**2)

def arc_length(x,ds,itermax=100,tol=1e-10):
    

    ds_old     = 0*ds
    grad_arc   = 0*ds
    for k in range(itermax):


        grad_arc_old = grad_arc
        grad_arc     = mvc_integrand_jacobian_arc_length(x,ds)
        
        #----------------------------------------------------#
        # Compute iteration step size:
        #----------------------------------------------------#
        dds = ds - ds_old
        dgrad_arc = grad_arc - grad_arc_old
        gamma_arc = 1e-10
        if dgrad_arc.dot(dgrad_arc)>0 and k>0:
            gamma_arc = np.abs(dds.dot(dgrad_arc))/dgrad_arc.dot(dgrad_arc)

        #----------------------------------------------------#
        # Update the solution:
        #----------------------------------------------------#
        ds_old = ds
        ds     = ds - gamma_arc*grad_arc

        ds_rel = np.sum(np.abs(grad_arc)/np.abs(ds))
        if ds_rel<tol:
            print(f"arc length iteration done, iter = {k:<4}, ls = {np.sum(ds):<10.10}")
            break


    return ds