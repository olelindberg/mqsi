from circle_from_three_points import circle_from_three_points

# Example:
p1 = (0, 0)
p2 = (1, 0)
p3 = (0, 1)

center, radius = circle_from_three_points(p1, p2, p3)
print(f"Center: {center}, Radius: {radius}")