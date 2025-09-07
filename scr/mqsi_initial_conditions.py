import numpy as np

from arc_parameter_circle import arc_parameter_circle
from circle_from_three_points import circle_from_three_points

from finite_difference_coefficients import finite_difference_coefficients




def mqsi_initial_conditions(curve,sol,center,radius,angle,points=[]):

    print(f"Setting equality constraints for {curve} ...")

    x,y,x_s,y_s,x_ss,y_ss = arc_parameter_circle(center,radius,angle)

    if curve=="wicket2":
        
        sol[2]  = x_ss[0]
        sol[5]  = y_ss[0]

        sol[8]  = x_ss[1]
        sol[11] = y_ss[1]

    elif curve=="wicket3":

        #  0:2  -  3:5
        #  6:8  -  9:11
        # 12:14 - 15:17

        # Vertex 1:
        sol[2]    = x_ss[0]
        sol[5]    = y_ss[0]

        # Vertex 2:
        sol[6+1]  = x_s [1]
        sol[6+2]  = x_ss[1]
        sol[6+4]  = y_s [1]
        sol[6+5]  = y_ss[1]

        # Vertex 3:
        sol[12+2] = x_ss[2]
        sol[12+5] = y_ss[2]

    elif curve=="wicket4":

        points                      = [[1.0, 0.5],[0.5, 1.0],[0.0, 0.5]]
        center,radius,arc_angles    = circle_from_three_points(points[0], points[1], points[2])
        x,y,x_s,y_s,x_ss,y_ss       = arc_parameter_circle(center,radius,arc_angles)

        # Vertex 1:
        sol[2]    = x_ss[0]
        sol[5]    = y_ss[0]

        # Vertex 2:
        sol[6+1]  = x_s [1]
        sol[6+2]  = x_ss[1]
        sol[6+4]  = y_s [1]
        sol[6+5]  = y_ss[1]

        # Vertex 3:
        sol[12+2] = x_ss[2]
        sol[12+5] = y_ss[2]

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
                sol[offset+2]  = x_ss[0]
                sol[offset+5]  = y_ss[0]
            elif i == len(points)-1: # Last point
                sol[offset+2] = x_ss[2]
                sol[offset+5] = y_ss[2]
            else:                   # Middle points
                sol[offset+1]  = x_s [1]
                sol[offset+2]  = x_ss[1]
                sol[offset+4]  = y_s [1]
                sol[offset+5]  = y_ss[1]

            offset += 6

    elif curve=="pointset":


        offset   = 0        
        for i in range(len(points)):
    
            
            if i == 0:               # First point


                x0 = points[i  ][0]
                x1 = points[i+1][0]
                x2 = points[i+2][0]

                y0 = points[i  ][1]
                y1 = points[i+1][1]
                y2 = points[i+2][1]

                dx0 = x1 - x0
                dx1 = x2 - x1

                dy0 = y1 - y0
                dy1 = y2 - y1

                ds0 = np.sqrt(dx0**2 + dy0**2) 
                ds1 = np.sqrt(dx1**2 + dy1**2)

                s0 = 0.0
                s1 = ds0
                s2 = ds0 + ds1 

                s = [s0,s1,s2]
                x = [x0,x1,x2]
                y = [y0,y1,y2]

                c2 = finite_difference_coefficients(s,s0,2)

                x_ss = c2.dot(x)
                y_ss = c2.dot(y)

                sol[offset+2] = x_ss
                sol[offset+5] = y_ss

            elif i == len(points)-1: # Last point

                x0 = points[i-2][0]
                x1 = points[i-1][0]
                x2 = points[i  ][0]

                y0 = points[i-2][1]
                y1 = points[i-1][1]
                y2 = points[i  ][1]

                dx0 = x1 - x0
                dx1 = x2 - x1

                dy0 = y1 - y0
                dy1 = y2 - y1

                ds0 = np.sqrt(dx0**2 + dy0**2) 
                ds1 = np.sqrt(dx1**2 + dy1**2)

                s0 = 0.0
                s1 = ds0
                s2 = ds0 + ds1 


                s = [s0,s1,s2]
                x = [x0,x1,x2]
                y = [y0,y1,y2]

                c2 = finite_difference_coefficients(s,s2,2)

                x_ss = c2.dot(x)
                y_ss = c2.dot(y)

                sol[offset+2] = x_ss
                sol[offset+5] = y_ss

            else:                   # Middle points

                x0 = points[i-1][0]
                x1 = points[i  ][0]
                x2 = points[i+1][0]

                y0 = points[i-1][1]
                y1 = points[i  ][1]
                y2 = points[i+1][1]

                dx0 = x1 - x0
                dx1 = x2 - x1

                dy0 = y1 - y0
                dy1 = y2 - y1

                ds0 = np.sqrt(dx0**2 + dy0**2) 
                ds1 = np.sqrt(dx1**2 + dy1**2)

                s0 = 0.0
                s1 = ds0
                s2 = ds0 + ds1 

                s = [s0,s1,s2]
                x = [x0,x1,x2]
                y = [y0,y1,y2]

                c1 = finite_difference_coefficients(s,s1,1)
                c2 = finite_difference_coefficients(s,s1,2)

                x_s  = c1.dot(x)
                x_ss = c2.dot(x)

                y_s  = c1.dot(y)
                y_ss = c2.dot(y)

                sol[offset+1]  = x_s
                sol[offset+2]  = x_ss
                sol[offset+4]  = y_s 
                sol[offset+5]  = y_ss

            offset += 6



    return sol