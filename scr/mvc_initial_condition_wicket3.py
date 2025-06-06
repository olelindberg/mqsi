import numpy as np
from circle_from_three_points import circle_from_three_points
from arc_parameter_circle import arc_parameter_circle
from constants import constants

def mvc_initial_condition_wicket3(r):

    print(f"Setting equality constraints for wicket3 ...")
    print("radius    = ",r)

    # Create the array of points:
    points     = np.array([[r, 0],[0,r],[-r,0]])
    num_points = len(points)

    center,radius = circle_from_three_points(points[0], points[1], points[2])

    # Calculate angles:
    angles = np.zeros(3)
    for j in range(3):
        dx = points[j] - center
        angles[j] = np.arctan2(dx[1], dx[0])

    # Arc parameter in meters [m]:
    s = radius*(angles-angles[0])

    # Initialize array for initial conditions:
    x0 = np.zeros((num_points,constants.NODE_DOFS))

    for i in range(num_points):

        # Calculate curve derivatives:
        x, y, x_s, y_s, x_ss, y_ss = arc_parameter_circle(radius,center,angles[0],s[-1],s[i])

        # Set initial condition:
        x0[i,0] = x
        x0[i,1] = x_s
        x0[i,2] = x_ss
        x0[i,3] = y
        x0[i,4] = y_s
        x0[i,5] = y_ss



    return x0,s,points