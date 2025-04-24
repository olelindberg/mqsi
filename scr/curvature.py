import numpy as np

def curvature(x_t,y_t,x_tt,y_tt):
    """
    Calculate the curvature of a curve defined by its parametric equations.
    The curvature is defined as the rate of change of the tangent vector
    with respect to the arc length. The formula used is:
    k = (x_t * y_tt - y_t * x_tt) / (x_t**2 + y_t**2)**(3/2)
    where:
    - x_t is the derivative of x with respect to t
    - y_t is the derivative of y with respect to t
    - x_tt is the second derivative of x with respect to t
    - y_tt is the second derivative of y with respect to t
    """
    eps = np.finfo(float).eps
    denom = (x_t**2 + y_t**2)**1.5
    denom[denom<eps] = eps  # Avoid division by zero
    return (x_t * y_tt - y_t * x_tt)/denom