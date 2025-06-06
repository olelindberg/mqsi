import numpy as np

def arc_parameter_circle(center,radius,theta):

    x    = center[0] + radius * np.cos(theta)
    y    = center[1] + radius * np.sin(theta)

    x_s  = -np.sin(theta)
    y_s  =  np.cos(theta)

    x_ss = -1/radius * np.cos(theta)
    y_ss = -1/radius * np.sin(theta)

    return x, y, x_s, y_s, x_ss, y_ss


