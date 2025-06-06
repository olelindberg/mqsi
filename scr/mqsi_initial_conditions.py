import numpy as np

from arc_parameter_circle import arc_parameter_circle

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



    return x0