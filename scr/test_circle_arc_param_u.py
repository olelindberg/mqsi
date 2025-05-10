import numpy as np
import matplotlib.pyplot as plt
from circle_arc_param_u import circle_arc_param_u

# Example parameters
r = 5
theta0 = np.pi / 6  # 30 degrees
delta_theta = np.pi / 2  # 90 degrees
theta_eval = np.pi/4

x, y, x_u, y_u, x_uu, y_uu = circle_arc_param_u(r, theta0, delta_theta,theta_eval)

plt.plot(x, y, label="Arc")
plt.quiver(x, y, x_u, y_u, color="red", scale=1, label="First derivative")
plt.quiver(x, y, x_uu, y_uu, color="blue", scale=1, label="Second derivative")
plt.axis('equal')
plt.legend()
plt.title("Circle Arc Parameterized by u ∈ [-1, 1]")
plt.show()