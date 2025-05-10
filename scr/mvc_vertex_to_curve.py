import numpy as np

def mvc_vertex_to_curve(x,i):

    cx      = np.zeros(6)

    vx    = x[i  ][0]       
    vx_s  = x[i  ][1]  
    vx_ss = x[i  ][2]  

    cx[0] = vx
    cx[1] = vx_s
    cx[2] = vx_ss
    
    vx    = x[i+1][0]    
    vx_s  = x[i+1][1]   
    vx_ss = x[i+1][2]  

    cx[3] = vx
    cx[4] = vx_s
    cx[5] = vx_ss

    cy = np.zeros(6)

    y    = x[i  ][3]
    y_s  = x[i  ][4]
    y_ss = x[i  ][5]

    cy[0] = y      
    cy[1] = y_s  
    cy[2] = y_ss
    
    y    = x[i+1][3]
    y_s  = x[i+1][4]
    y_ss = x[i+1][5]

    cy[3] = y    
    cy[4] = y_s  
    cy[5] = y_ss

    return cx,cy
