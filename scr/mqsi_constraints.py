import numpy as np
from arc_parameter_circle import arc_parameter_circle


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


    return bc_dof,bc_value