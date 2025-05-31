import numpy as np

def circle_from_three_points(p1, p2, p3):
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

    if abs(a) < 1e-10:
        raise ValueError("Points are colinear")

    Dx = np.linalg.det(np.hstack([B, A[:, [1, 2]]]))
    Dy = np.linalg.det(np.hstack([A[:, [0]], B, A[:, [2]]]))
    C  = np.linalg.det(np.hstack([A[:, [0, 1]], B]))

    cx = 0.5 * Dx / a
    cy = 0.5 * Dy / a
    r  = np.sqrt(cx**2 + cy**2 + C / a)

    center = np.array([cx, cy])
    radius = r
    
    return center, radius


