from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt

from mvc_integrand_jacobian        import mvc_integrand_jacobian
from mvc_objective_function        import mvc_objective_function
from mvc_initial_condition_wicket3 import mvc_initial_condition_wicket3
from mvc_vertex_to_curve           import mvc_vertex_to_curve
from hermite_quintic               import hermite_quintic
from curvature                     import curvature
from circle_from_three_points      import circle_from_three_points
from circle_arc_param_u            import circle_arc_param_u
from constants                     import constants

# node dof arragement
# 0 - x
# 1 - x'
# 2 - x''
# 3 - y
# 4 - y'
# 5 - y''




def assign_constraints(x0, constraint_ids, constraint_values):
    for i in range(len(constraint_ids)):
        for j in range(len(constraint_ids[i])):
            k = i*constants.NODE_DOFS + constraint_ids[i][j]
            x0[k] = constraint_values[i][j]
    return x0

def assign_constraints_grad(x0, constraint_ids):
    for i in range(len(constraint_ids)):
        for j in range(len(constraint_ids[i])):
            k = i*constants.NODE_DOFS + constraint_ids[i][j]
            x0[k] = 0.0
    return x0


maxiter     = 10
tol         = 1e-5
solver_type = "gradient_descent" # "trust-constr" # "SLSQP" # "L-BFGS-B" # "dogleg" # "trust-ncg"
curve       = "wicket3"
show_figures = True


if curve == "wicket2":

    # Curve parameters:
    r =1.0

    print(f"Setting equality constraints for {curve} ...")
    print("radius    = ",r)

    # Create the array of points:
    points     = np.array([[r, 0],[-r,0]])
    num_points = len(points)

    # Create the circle:
    center,radius = circle_from_three_points(points[0], np.array([0,r]), points[1])

    # Calculate angles:
    angles = np.zeros(num_points)
    for i in range(num_points):
        dx = points[i] - center
        angles[i] = np.arctan2(dx[1], dx[0])

    # Calculate arc length of curve:
    arc_length = angles[1] - angles[0]
    if arc_length < 0:
        arc_length += 2*np.pi
    
    # Initialize array for constraint values:
    constraint_values = np.zeros((num_points,4))

    # Initialize array for initial conditions:
    x0 = np.zeros((num_points,constants.NODE_DOFS))

    for i in range(num_points):

        # First point:
        if i==0:

            # Calculate curve derivatives:
            x, y, x_u, y_u, x_uu, y_uu = circle_arc_param_u(radius,center,angles[0],arc_length,angles[i])

            # Add constraints:
            constraint_values[i,0] = x
            constraint_values[i,1] = x_u
            constraint_values[i,2] = y
            constraint_values[i,3] = y_u

        elif i==num_points-1:

            # Calculate curve derivatives:
            x, y, x_u, y_u, x_uu, y_uu = circle_arc_param_u(radius,center,angles[0],arc_length,angles[i])

            # Add constraints:
            constraint_values[i,0] = x
            constraint_values[i,1] = x_u
            constraint_values[i,2] = y
            constraint_values[i,3] = y_u
        
        else:
            NotImplementedError("Only implemented for boundary points")


        # Set initial condition:
        x0[i,0] = x
        x0[i,1] = x_u
        x0[i,2] = x_uu
        x0[i,3] = y
        x0[i,4] = y_u
        x0[i,5] = y_uu


    constraint_ids = [[0, 1, 3, 4],[ 0, 1, 3, 4]]

    x0 = x0.flatten()    

elif curve == "wicket3":

    radius      = 1
    x0,s,points = mvc_initial_condition_wicket3(radius)
    num_points  = len(points)

    # Initialize array for constraint values:
    constraint_values = []
    constraint_ids    = []

    for i in range(num_points):

        x   = x0[i,:][0]
        x_u = x0[i,:][1]
        y   = x0[i,:][3] 
        y_u = x0[i,:][4] 

        if 0<i and i<num_points-1:
            # Add constraints:
            cv = [x,y]
            ci = [0,3]
        else:
            # Add constraints:
            cv = [x,x_u,y,y_u]
            ci = [0,1,3,4]

        # Append constraint values:
        constraint_values.append(cv)
        constraint_ids.append(ci)

    x0 = x0.flatten()

elif curve == "points2":
    constraint_ids    = [[  0,   1,   3,   4],[  0,   1,   3,   4]]
    constraint_values  = [[0.0, 1.0, 0.0, 0.0],[1.0, 1.0, 1.0, 0.0]]
elif curve == "points3":
    constraint_ids    = [[  0,   3],[0,     3],[0, 3]]
    constraint_values  = [[0.0, 0.0],[0.5,0.25],[1, 1]]


#-----------------------------------------------------------------------------#
# Equality constraints
#-----------------------------------------------------------------------------#
equality_constraints = []
for i in range(len(constraint_ids)):
    for j in range(len(constraint_ids[i])):
        k = i*constants.NODE_DOFS + constraint_ids[i][j]
        eq = {'type': 'eq', 'fun': lambda x,i=i,j=j,k=k: x[k] - constraint_values[i][j]}
        equality_constraints.append(eq)

#-----------------------------------------------------------------------------#
# Initial condition
#-----------------------------------------------------------------------------#
x0 = assign_constraints(x0, constraint_ids, constraint_values)

if curve=="points3":

    x0  = [[  0, 0.5,  1,    0, 0.25,  1],
           [0.5,   1,  1, 0.25,    1,  1],
           [  1, 0.5,  1,    1, 0.75,  1]]
    x0 = np.array(x0).flatten()

    #x0[7]  = 1 # derivative at the middle point
    #x0[10] = 1 # derivative at the middle point


#-----------------------------------------------------------------------------#
# Solve:
#-----------------------------------------------------------------------------#
if solver_type == "trust-constr":

    print("Solving with trust-constr ...")
    
    options = {"maxiter" : maxiter,'disp': True, "verbose" : 1}
    method =  'trust-constr' # 'SLSQP' #
    
    res     = minimize(mvc_objective_function, x0, jac=mvc_integrand_jacobian,constraints=equality_constraints , tol=tol, method=method,options=options)

    x = res.x

elif solver_type == "gradient_descent":

    print("Solving with gradient descent ...")

    #----------------------------------------------------#
    # Initialization of variables:
    #----------------------------------------------------#
    x     = x0
    x_old = 0*x0
    grad  = np.zeros(len(x))

    #----------------------------------------------------#
    # Initialization of norms:
    #----------------------------------------------------#
    x_norm_init = np.linalg.norm(x)
    x_norm_old  = x_norm_init
    
    #----------------------------------------------------#
    # Main loop:
    #----------------------------------------------------#
    for i in range(maxiter):

        #----------------------------------------------------#
        # Compute the gradient:
        #----------------------------------------------------#
        grad_old = grad
        grad     = mvc_integrand_jacobian(x,s)
        grad     = assign_constraints_grad(grad, constraint_ids)

        #----------------------------------------------------#
        # Compute iteration step size:
        #----------------------------------------------------#
        dx    = x - x_old
        dgrad = grad - grad_old
        gamma = 1e-10
        if dgrad.dot(dgrad)>0:
            gamma = np.abs(dx.dot(dgrad))/dgrad.dot(dgrad)

        #----------------------------------------------------#
        # Update the solution:
        #----------------------------------------------------#
        x_old = x
        x     = x - gamma*grad
        
        #----------------------------------------------------#
        # Check convergence:
        #----------------------------------------------------#
        x_norm     = np.linalg.norm(x)
        dx_norm    = np.abs(x_norm - x_norm_old)/x_norm_init
        x_norm_old = x_norm
        print(f"Iteration {i:<3}  gamma: {gamma:<10.4e}  dx_norm: {dx_norm:<10.4e}")
        if dx_norm < tol:
            print(f"Converged in {i} iterations.")
            break

if show_figures:

    #------------------------------------------------------------#
    # Plotting:
    #------------------------------------------------------------#
    xi = np.arange(-1, 1.01, 0.01)

    x = np.reshape(x,(len(constraint_ids),6))

    for i in range(x.shape[0]-1):

        cx,cy = mvc_vertex_to_curve(x,i)


        s_xi = 0.5*(s[i+1] - s[i])
        H0  = hermite_quintic(0,xi,s_xi)
        Ht  = hermite_quintic(1,xi,s_xi)
        Htt = hermite_quintic(2,xi,s_xi)

        xx  = H0@cx
        yy  = H0@cy

        cx_t  = Ht@cx
        cy_t  = Ht@cy

        cx_tt  = Htt@cx
        cy_tt  = Htt@cy

        k  = curvature(cx_t, cy_t, cx_tt, cy_tt)

        plt.figure(1)
        plt.plot(xi, xx)
        plt.title('x')
        plt.grid(True)

        plt.figure(2)
        plt.plot(xi, yy)
        plt.title('y')
        plt.grid(True)

        plt.figure(3)
        plt.plot(xi, k)
        plt.title('curvature')
        plt.grid(True)

        plt.figure(4)
#        angle = np.arange(angles[0],angles[2],arc_length/100)
#        xxx,yyy,tmp,tmp,tmp,tmp = circle_arc_param_u(radius,center,angles[0],arc_length,angle)
#        plt.plot(xxx,yyy,'k')


        for i in range(num_points):

            plt.plot(points[i,0],points[i,1],'bo')

        plt.plot(xx, yy)
#        plt.plot(center[0],center[1],'ro')
        plt.title('x-y')
        plt.grid(True)

    plt.show()