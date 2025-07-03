import numpy as np

from arc_parameter_circle import arc_parameter_circle
from circle_from_three_points import circle_from_three_points

def mqsi_initial_conditions(curve,x0,center,radius,angle):

    print(f"Setting equality constraints for {curve} ...")

    x,y,x_s,y_s,x_ss,y_ss = arc_parameter_circle(center,radius,angle)

    if curve=="wicket2":
        
        x0[2]  = x_ss[0]
        x0[5]  = y_ss[0]

        x0[8]  = x_ss[1]
        x0[11] = y_ss[1]

    elif curve=="wicket3":

        #  0:2  -  3:5
        #  6:8  -  9:11
        # 12:14 - 15:17

        # Vertex 1:
        x0[2]    = x_ss[0]
        x0[5]    = y_ss[0]

        # Vertex 2:
        x0[6+1]  = x_s [1]
        x0[6+2]  = x_ss[1]
        x0[6+4]  = y_s [1]
        x0[6+5]  = y_ss[1]

        # Vertex 3:
        x0[12+2] = x_ss[2]
        x0[12+5] = y_ss[2]

    elif curve=="wicket4":

        points                      = [[1.0, 0.5],[0.5, 1.0],[0.0, 0.5]]
        center,radius,arc_angles    = circle_from_three_points(points[0], points[1], points[2])
        x,y,x_s,y_s,x_ss,y_ss       = arc_parameter_circle(center,radius,arc_angles)

        # Vertex 1:
        x0[2]    = x_ss[0]
        x0[5]    = y_ss[0]

        # Vertex 2:
        x0[6+1]  = x_s [1]
        x0[6+2]  = x_ss[1]
        x0[6+4]  = y_s [1]
        x0[6+5]  = y_ss[1]

        # Vertex 3:
        x0[12+2] = x_ss[2]
        x0[12+5] = y_ss[2]

    elif curve=="wicket5":

        points = [[1.0, 0.0],[0.0, 1.0],[-1.0,0.0 ],[0.0,-1.0]]

        offset   = 0        
        for i in range(len(points)):
    
            if i == 0:               # First point
                center, radius,arc_angles = circle_from_three_points(points[i], points[i+1], points[i+2])
            elif i == len(points)-1: # Last point
                center, radius,arc_angles = circle_from_three_points(points[i-2], points[i-1], points[i])   
            else:                    # Middle points
                center, radius,arc_angles = circle_from_three_points(points[i-1], points[i], points[i+1])

            x,y,x_s,y_s,x_ss,y_ss = arc_parameter_circle(center,radius,arc_angles)

            if i == 0:               # First point
                x0[offset+2]  = x_ss[0]
                x0[offset+5]  = y_ss[0]
            elif i == len(points)-1: # Last point
                x0[offset+2] = x_ss[2]
                x0[offset+5] = y_ss[2]
            else:                   # Middle points
                x0[offset+1]  = x_s [1]
                x0[offset+2]  = x_ss[1]
                x0[offset+4]  = y_s [1]
                x0[offset+5]  = y_ss[1]

            offset += 6


    return x0