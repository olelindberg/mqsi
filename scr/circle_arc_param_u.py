import numpy as np
import matplotlib.pyplot as plt

def circle_arc_param_u(radius,center,theta0,arc_length,s):

    theta = s/radius + theta0

    x = center[0] + radius * np.cos(theta)
    y = center[1] + radius * np.sin(theta)

    x_u = -np.sin(theta)
    y_u =  np.cos(theta)

    x_uu = -1/radius * np.cos(theta)
    y_uu = -1/radius * np.sin(theta)

    return x, y, x_u, y_u, x_uu, y_uu


