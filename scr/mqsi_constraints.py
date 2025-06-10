import numpy as np

from arc_parameter_circle     import arc_parameter_circle
from circle_from_three_points import circle_from_three_points

def mqsi_constraints(curve,center,radius,theta):

    print(f"Setting equality constraints for {curve} ...")

    x,y,x_s,y_s,x_ss,y_ss = arc_parameter_circle(center,radius,theta)

    if curve == "wicket2":
                     
        bc_dof    = [[   0,      1,    3,      4],[   0,      1,    3,      4]]
        bc_value  = [[x[0], x_s[0], y[0], y_s[0]],[x[1], x_s[1], y[1], y_s[1]]]

    elif curve == "wicket3":

        bc_dof   = []
        bc_value = []

        bc_dof.append  ([   0,      1,    3,      4])
        bc_value.append([x[0], x_s[0], y[0], y_s[0]])

        bc_dof.append  ([   0,    3])
        bc_value.append([x[1], y[1]])

        bc_dof.append  ([   0,      1,    3,      4])
        bc_value.append([x[2], x_s[2], y[2], y_s[2]])


    elif curve == "wicket4":

        p1 = (1.0, 0.5) 
        p2 = (0.5, 1.0)
        p3 = (0.0, 0.5)

        center, radius,arc_angles = circle_from_three_points(p1, p2, p3)
        x,y,x_s,y_s,x_ss,y_ss = arc_parameter_circle(center,radius,arc_angles)

        bc_dof   = []
        bc_value = []

        bc_dof.append  ([   0,      1,    3,      4])
        bc_value.append([x[0], x_s[0], y[0], y_s[0]])

        bc_dof.append  ([   0,    3])
        bc_value.append([x[1], y[1]])

        bc_dof.append  ([   0,      1,    3,      4])
        bc_value.append([x[2], x_s[2], y[2], y_s[2]])

    elif curve == "wicket4":

        points = [[1.0, 0.5],[0.5, 1.0],[0.0, 0.5]]
        center, radius,arc_angles = circle_from_three_points(points[0], points[1], points[2])
        x,y,x_s,y_s,x_ss,y_ss = arc_parameter_circle(center,radius,arc_angles)

        bc_dof   = []
        bc_value = []

        bc_dof.append  ([   0,      1,    3,      4])
        bc_value.append([x[0], x_s[0], y[0], y_s[0]])

        bc_dof.append  ([   0,    3])
        bc_value.append([x[1], y[1]])

        bc_dof.append  ([   0,      1,    3,      4])
        bc_value.append([x[2], x_s[2], y[2], y_s[2]])

    elif curve == "wicket5":

        points = [[1.0, 0.5],[0.5, 1.0],[0.0, 0.5],[0.5, 0.0]]

        bc_dof   = []
        bc_value = []

        for i in range(len(points)):
            print(f"Processing point {i} of {len(points)} ...")
            if i == 0:               # First point
                center, radius,arc_angles = circle_from_three_points(points[i], points[i+1], points[i+2])
            elif i == len(points)-1: # Last point
                center, radius,arc_angles = circle_from_three_points(points[i-2], points[i-1], points[i])   
            else:                    # Middle points
                center, radius,arc_angles = circle_from_three_points(points[i-1], points[i], points[i+1])

            x,y,x_s,y_s,x_ss,y_ss = arc_parameter_circle(center,radius,arc_angles)

            if i == 0:               # First point
                bc_dof.append  ([   0,      1,    3,      4])
                bc_value.append([x[0], x_s[0], y[0], y_s[0]])
            elif i == len(points)-1: # Last point
                bc_dof.append  ([   0,      1,    3,      4])
                bc_value.append([x[2], x_s[2], y[2], y_s[2]])
            else:                    # Middle points
                bc_dof.append  ([   0,    3])
                bc_value.append([x[1], y[1]])

    return bc_dof,bc_value