import numpy as np

from circle_from_three_points import circle_from_three_points

# Example:
p1 = (1.0, 0.5)
p2 = (0.5, 1.0)
p3 = (0.0, 0.5)

center, radius,arc_angles = circle_from_three_points(p1, p2, p3)

print(f"Center: ({center[0]}, {center[1]}), Radius: {radius},  arc_angles: {arc_angles}")



