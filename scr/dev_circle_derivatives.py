import numpy as np
import matplotlib.pyplot as plt

def circle_arc_param_u(r, theta0, delta_theta, num_points=100):

    L     = r * delta_theta
    u     = np.linspace(-1, 1, num_points)
    s     = (L / 2) * (u + 1)
    theta = s / r + theta0

    x     = r * np.cos(theta)
    y     = r * np.sin(theta)
    x_u   = (L / 2) * (-np.sin(theta))
    y_u   = (L / 2) * ( np.cos(theta))
    x_uu = (L/2)**2 * (-1/r) * np.cos(theta)
    y_uu = (L/2)**2 * (-1/r) * np.sin(theta)

    return x, y, x_u, y_u, x_uu, y_uu

# Example parameters
r = 5
theta0 = np.pi / 6  # 30 degrees
delta_theta = np.pi / 2  # 90 degrees

x, y, x_u, y_u, x_uu, y_uu = circle_arc_param_u(r, theta0, delta_theta)
print(x_uu, y_uu)
plt.plot(x, y, label="Arc")
plt.quiver(x, y, x_u, y_u, color="red", scale=1, label="First derivative")
plt.quiver(x, y, x_uu, y_uu, color="blue", scale=1, label="Second derivative")
plt.axis('equal')
plt.legend()
plt.title("Circle Arc Parameterized by u ∈ [-1, 1]")
plt.show()
