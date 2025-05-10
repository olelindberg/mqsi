import numpy as np
import matplotlib.pyplot as plt

def circle_arc_param_u(radius,center,theta0,arc_length,s):
    """
    Parametrizes a circular arc using the parameter u ∈ [-1, 1].
    
    Parameters:
    r : float
        Radius of the circle.
    theta0 : float
        Initial angle in radians.
    delta_theta : float
        Angle subtended by the arc in radians.
    num_points : int
        Number of points to generate along the arc.
    
    Returns:
    x : np.ndarray
        x-coordinates of the points on the arc.
    y : np.ndarray
        y-coordinates of the points on the arc.
    x_u : np.ndarray
        First derivative of x with respect to u.
    y_u : np.ndarray
        First derivative of y with respect to u.
    x_uu : np.ndarray
        Second derivative of x with respect to u.
    y_uu : np.ndarray
        Second derivative of y with respect to u.
    
    """
    theta = s/radius + theta0

    x = center[0] + radius * np.cos(theta)
    y = center[1] + radius * np.sin(theta)

    x_u = -np.sin(theta)
    y_u =  np.cos(theta)

    x_uu = -1/radius * np.cos(theta)
    y_uu = -1/radius * np.sin(theta)

    return x, y, x_u, y_u, x_uu, y_uu


