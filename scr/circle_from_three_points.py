import numpy as np

def arc_parameter(center,point):

    dx = point[0] - center[0]
    dy = point[1] - center[1]
    theta = np.arctan2(dy, dx)
    return theta

def circle_from_three_points(p1, p2, p3):

    scale = np.max([np.max(p1),np.max(p2),np.max(p3)])

    p1 = 1/scale*p1
    p2 = 1/scale*p2
    p3 = 1/scale*p3


    A = np.array([
        [p1[0], p1[1], 1],
        [p2[0], p2[1], 1],
        [p3[0], p3[1], 1]
    ])
    
    B = np.array([
        [p1[0]**2 + p1[1]**2],
        [p2[0]**2 + p2[1]**2],
        [p3[0]**2 + p3[1]**2]
    ])

    # Determinants
    a = np.linalg.det(A)
    center = np.zeros(2)
    radius = 0

    #print(a)

    if abs(a) < 1e-10:
        return p2,np.finfo(p2.dtype).max,0

    Dx = np.linalg.det(np.hstack([B, A[:, [1, 2]]]))
    Dy = np.linalg.det(np.hstack([A[:, [0]], B, A[:, [2]]]))
    C  = np.linalg.det(np.hstack([A[:, [0, 1]], B]))

    cx = 0.5 * Dx / a
    cy = 0.5 * Dy / a
    r  = np.sqrt(cx**2 + cy**2 + C / a)

    center = scale*np.array([cx, cy])
    radius = scale*r

    # Calculate the angle of the center point relative to the first point
    arc_angles = np.array([arc_parameter(center, p1),arc_parameter(center, p2),arc_parameter(center, p3)])

    return center,radius,arc_angles


