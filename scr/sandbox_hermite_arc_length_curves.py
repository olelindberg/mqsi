from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

from mvc_integrand_jacobian     import mvc_integrand_jacobian
from mvc_objective_function     import mvc_objective_function
from hermite_quintic            import hermite_quintic
from hermite_quintic_derivatives import hermite_quintic_derivatives
from curvature                  import curvature
from circle_from_three_points   import circle_from_three_points
from circle_arc_param_u         import circle_arc_param_u
from mvc_initial_condition_wicket3 import mvc_initial_condition_wicket3
from mvc_vertex_to_curve      import mvc_vertex_to_curve

# node dof arragement
# 0 - x
# 1 - x'
# 2 - x''
# 3 - y
# 4 - y'
# 5 - y''







@dataclass(frozen=True)
class Constants:
    NODE_DOFS: int = 6

curve       = "wicket3"
#curve       = "points3"
show_figures = True

if curve == "wicket3":

    radius      = 5
    x0,s,points = mvc_initial_condition_wicket3(radius)
    num_points  = len(points)

elif curve == "points3":
    num_points = 3
    points  = np.array([[  0,    0, ],
                        [0.5, 0.25, ],
                        [  1,    1, ]])
    x0  = [[  0, 0.5,  1,    0, 0.25,  1],
           [0.5,   1,  1, 0.25,    1,  1],
           [  1, 0.5,  1,    1, 0.75,  1]]
    x0 = np.array(x0).flatten()

#-----------------------------------------------------------------------------#
# Initial condition
#-----------------------------------------------------------------------------#
x = x0

if show_figures:

    #------------------------------------------------------------#
    # Plotting:
    #------------------------------------------------------------#
    t     = np.arange(0, 1.1, 0.1)

    x = np.reshape(x,(num_points,6))
    
    fig, ax = plt.subplots(2, 2, figsize=(10, 4))

    for i in range(x.shape[0]-1):

        cx,cy = mvc_vertex_to_curve(x,i)

        s_t = (s[i+1] - s[i])        
        xx, yy, cx_s, cy_s, cx_ss, cy_ss = hermite_quintic_derivatives(s_t,t,cx,cy)

        # Curvature:
        k = np.sqrt(cx_ss**2 + cy_ss**2)

        ax[0,0].plot(t, cx_s,label="cx_s " + str(i))
        ax[0,0].plot(t, cy_s,label="cy_s " + str(i))
        ax[0,0].set_title('c_s')
        ax[0,0].grid(True)
        ax[0,0].legend()

        ax[1,0].plot(t, cx_ss,label="cx_ss " + str(i))
        ax[1,0].plot(t, cy_ss,label="cy_ss " + str(i))
        ax[1,0].set_title('c_ss')
        ax[1,0].grid(True)
        ax[1,0].legend()

        ax[0,1].plot(t, k)
        ax[0,1].set_title('curvature')
        ax[0,1].grid(True)
        
        ax[1,1].plot(points[:,0],points[:,1],'bo')
        ax[1,1].plot(xx, yy,label=str(i))
        ax[1,1].axis('equal')
        ax[1,1].set_title('x-y')
        ax[1,1].grid(True)
        ax[1,1].legend()

    plt.show()