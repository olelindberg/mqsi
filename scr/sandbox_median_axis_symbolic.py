import sympy as sp

# Define the parametric equations for the boundary curve C(t)
t = sp.symbols('t')
x_t = sp.cos(t)  # Example: Parametric curve (circle)
y_t = sp.sin(t)

# Calculate the tangent vector
x_t_derivative = sp.diff(x_t, t)
y_t_derivative = sp.diff(y_t, t)

# Calculate the normal vector (rotate tangent vector by 90 degrees)
normal_x = -y_t_derivative
normal_y = x_t_derivative

# Parametrize the normal at two points (t1, t2)
t1, t2 = sp.symbols('t1 t2')
p1_x = x_t.subs(t, t1) + normal_x.subs(t, t1)
p1_y = y_t.subs(t, t1) + normal_y.subs(t, t1)

p2_x = x_t.subs(t, t2) + normal_x.subs(t, t2)
p2_y = y_t.subs(t, t2) + normal_y.subs(t, t2)

# Solve for intersection point of normals
equations = [p1_x - p2_x, p1_y - p2_y]
print("Equations to solve for intersection points:")
for eq in equations:
    print(eq)
intersection_points = sp.solve(equations, (t1, t2))

# Now `intersection_points` gives us exact locations for the centers of maximal balls
print(intersection_points)


a1 = sp.symbols('a')
a2 = sp.symbols('a')
a3 = sp.symbols('a')

b1 = sp.symbols('b')
b2 = sp.symbols('b')
b3 = sp.symbols('b')

t1 = sp.symbols('t1')
t2 = sp.symbols('t2')

p = sp.symbols('p')

p1 = a1 + a2*t1 + a3*t1**2
p2 = b1 + b2*t2 + b3*t2**2

v1 = p - p1
v2 = p - p2

f1 = sp.sqrt(v1**2 + v1[1]**2)



p