import numpy as np
import matplotlib.pyplot as plt

def arc_parameter_circle(r, theta0, dtheta, num_points=100):

    u     = np.linspace(0, 1, num_points)
    theta = dtheta*u + theta0

    x    =            r*np.cos(theta)
    y    =            r*np.sin(theta)
    x_u  =    -r*dtheta*np.sin(theta)
    y_u  =     r*dtheta*np.cos(theta)
    x_uu = -r*dtheta**2*np.cos(theta)
    y_uu = -r*dtheta**2*np.sin(theta)

    return x, y, x_u, y_u, x_uu, y_uu

# Example parameters
r = 5
theta0 = np.pi / 6  # 30 degrees
delta_theta = np.pi / 2  # 90 degrees

x, y, x_u, y_u, x_uu, y_uu = arc_parameter_circle(r, theta0, delta_theta)
print(x_uu, y_uu)
plt.plot(x, y, label="Arc")
plt.quiver(x, y, x_u, y_u, color="red", scale=1, label="First derivative")
plt.quiver(x, y, x_uu, y_uu, color="blue", scale=1, label="Second derivative")
plt.axis('equal')
plt.legend()
plt.title("Circle Arc Parameterized by u ∈ [-1, 1]")
plt.show()
