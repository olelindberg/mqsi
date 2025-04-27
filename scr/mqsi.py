from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

from mvc_integrand_jacobian import mvc_integrand_jacobian
from mvc_objective_function import mvc_objective_function
from hermite_quintic        import hermite_quintic
from curvature              import curvature

# node dof arragement
# 0 - x
# 1 - x'
# 2 - x''
# 3 - y
# 4 - y'
# 5 - y''


@dataclass(frozen=True)
class Constants:
    NODE_DOFS: int = 6

def assign_constraints(x0, bc_dof, bc_value):
    for i in range(len(bc_dof)):
        for j in range(len(bc_dof[i])):
            k = i*Constants.NODE_DOFS + bc_dof[i][j]
            x0[k] = bc_value[i][j]
    return x0

def assign_constraints_grad(x0, bc_dof):
    for i in range(len(bc_dof)):
        for j in range(len(bc_dof[i])):
            k = i*Constants.NODE_DOFS + bc_dof[i][j]
            x0[k] = 0.0
    return x0


maxiter     = 1000
tol         = 1e-5
solver_type = "gradient_descent" # "trust-constr" # "SLSQP" # "L-BFGS-B" # "dogleg" # "trust-ncg"
curve       = "wicket3"
show_figures = False

if curve == "wicket2":

    r = 1.0
    l = np.pi*r
    k = 1/r

    print("Setting equality constraints for wicket2 ...")
    print("curve     = ",curve)
    print("radius    = ",r)
    print("length    = ",l)
    print("curvature = ",k)

    bc_dof    = [[0, 1, 3,         4],[ 0, 1, 3,          4]]
    bc_value  = [[r, 0, 0, r*np.pi/2],[-r, 0, 0, -r*np.pi/2]]

elif curve == "wicket3":

    r     = 10
    angle = np.pi/4
    k     = 1/r

#    cx    = [r,       0, -r*angle**2, 0, -r*angle,           0]
#    cy    = [0, r*angle,           0, 1,        0, -r*angle**2]

    bc_dof    = [[0,   1,   3,       4],[0,3],[ 0,   1,   3,        4]]
    bc_value  = [[r, 0.0, 0.0, r*angle],[0,r],[-r, 0.0, 0.0, -r*angle]]

elif curve == "points2":
    bc_dof    = [[  0,   1,   3,   4],[  0,   1,   3,   4]]
    bc_value  = [[0.0, 1.0, 0.0, 0.0],[1.0, 1.0, 1.0, 0.0]]
elif curve == "points3":
    bc_dof    = [[  0,   3],[0,     3],[0, 3]]
    bc_value  = [[0.0, 0.0],[0.5,0.25],[1, 1]]


#-----------------------------------------------------------------------------#
# Equality constraints
#-----------------------------------------------------------------------------#
equality_constraints = []
for i in range(len(bc_dof)):
    for j in range(len(bc_dof[i])):
        k = i*Constants.NODE_DOFS + bc_dof[i][j]
        eq = {'type': 'eq', 'fun': lambda x,i=i,j=j,k=k: x[k] - bc_value[i][j]}
        equality_constraints.append(eq)

#-----------------------------------------------------------------------------#
# Initial condition
#-----------------------------------------------------------------------------#
x0 = np.zeros(Constants.NODE_DOFS * len(bc_dof))
x0 = assign_constraints(x0, bc_dof, bc_value)

if curve=="wicket2":
    
    print("Setting initial condition for wicket2 ...")
    
    x0[2] = -r/4*np.pi**2
    x0[8] =  r/4*np.pi**2

elif curve=="wicket3":
    # 0     - 5
    # 6:8   - 9:11
    # 12:14 - 15:17

    # Vertex 1:
    x0[2]  = -r*angle**2

    # Vertex 2:
    x0[7]  = -r*angle
    x0[11] = -r*angle**2

    # Vertex 3:
    x0[14] =  r*angle**2

elif curve=="points3":

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
        grad     = mvc_integrand_jacobian(x)
        grad     = assign_constraints_grad(grad, bc_dof)

#        print("grad")
#        print(grad)

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
    H0  = hermite_quintic(0,xi)
    Ht  = hermite_quintic(1,xi)
    Htt = hermite_quintic(2,xi)

    x = np.reshape(x,(len(bc_dof),6))
    #print("x")
    #print(x)

    for i in range(x.shape[0]-1):

        cx = np.zeros(6)
        cx[0:3] = x[i][0:3]
        cx[3:6] = x[i+1][0:3]

        cy = np.zeros(6)
        cy[0:3] = x[i][3:6]
        cy[3:6] = x[i+1][3:6]

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
        plt.plot(xx, yy)
        plt.title('x-y')
        plt.grid(True)

    plt.show()