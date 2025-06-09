import numpy as np
from scipy.optimize import minimize
from matplotlib import pyplot as plt

# Parametric curve (circle for simplicity)
def parametric_curve(t):
    return np.array([2*np.cos(t), np.sin(t)])

# Normal vector at a given point
def normal_vector(t):
    return np.array([-2*np.sin(t), np.cos(t)])

# Distance function for optimization
def distance_to_boundary(p, t1, t2):
    p1 = parametric_curve(t1)
    p2 = parametric_curve(t2)
    # Compute distance from point to boundary points (for simplicity)
    v1 = p-p1
    v2 = p-p2
    d  = (np.linalg.norm(v2)-np.linalg.norm(v1))**2
    return  d

# Objective function to minimize: the distance to the curve
def objective(p, t1, t2):
    return distance_to_boundary(p, t1, t2)

# Constraints: the normals should point toward the point p
def constraint1(p, t1):
    p1 = parametric_curve(t1)
    normal = normal_vector(t1)
    return np.dot(normal, p - p1)

# Define initial guess for point p
initial_guess = np.array([0.1, 0.1])

# Set optimization bounds and constraints
t1, t2 = 0.0, np.pi/2  # For simplicity, pick two points on the curve
constraints = [{'type': 'ineq', 'fun': constraint1, 'args': (0,)},
               {'type': 'ineq', 'fun': constraint1, 'args': (np.pi/2,)},
               {'type': 'ineq', 'fun': constraint1, 'args': (np.pi  ,)},
               {'type': 'ineq', 'fun': constraint1, 'args': (2*np.pi,)}]

# Solve the optimization problem
result = minimize(objective, initial_guess, args=(t1, t2), constraints=constraints)
optimal_point = result.x
print("Optimal Point:", optimal_point)


angles = np.linspace(0, 2*np.pi, 100)
curve_points = parametric_curve(angles)
plt.plot(curve_points[0], curve_points[1], label='Curve')
plt.scatter(optimal_point[0], optimal_point[1], color='red', label='Optimal Point')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Optimal Point on Parametric Curve')
plt.legend()
plt.axis('equal')
plt.grid()
plt.show()
