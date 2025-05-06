import numpy as np

def circle_from_three_points_np(p1, p2, p3):
    
    A = np.array([
        [p1[0], p1[1], 1],
        [p2[0], p2[1], 1],
        [p3[0], p3[1], 1]
    ])
    
    D = np.linalg.det(A)
    if abs(D) < 1e-10:
        raise ValueError("Points are colinear")

    # Calculate terms for center
    a = np.array([
        [p1[0]**2 + p1[1]**2, p1[1], 1],
        [p2[0]**2 + p2[1]**2, p2[1], 1],
        [p3[0]**2 + p3[1]**2, p3[1], 1]
    ])
    d = np.linalg.det(a)

    b = np.array([
        [p1[0]**2 + p1[1]**2, p1[0], 1],
        [p2[0]**2 + p2[1]**2, p2[0], 1],
        [p3[0]**2 + p3[1]**2, p3[0], 1]
    ])
    e = np.linalg.det(b)

    c = np.array([
        [p1[0]**2 + p1[1]**2, p1[0], p1[1]],
        [p2[0]**2 + p2[1]**2, p2[0], p2[1]],
        [p3[0]**2 + p3[1]**2, p3[0], p3[1]]
    ])
    f = np.linalg.det(c)

    # Circle center (h, k) and radius r
    h =  0.5 * d / D
    k = -0.5 * e / D
    r = np.sqrt(h**2 + k**2 + f / D)

    return (h, k, r)

# Example:
p1 = (0, 0)
p2 = (1, 0)
p3 = (0, 1)

center_x, center_y, radius = circle_from_three_points_np(p1, p2, p3)
print(f"Center: ({center_x}, {center_y}), Radius: {radius}")