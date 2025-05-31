from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

from mvc_integrand_jacobian import mvc_integrand_jacobian
from mvc_objective_function import mvc_objective_function
from hermite_quintic        import hermite_quintic
from curvature              import curvature
from curve_metrics          import arc_length


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


maxiter     = 1
tol         = 1e-8
solver_type = "gradient_descent" # "trust-constr" # "SLSQP" # "L-BFGS-B" # "dogleg" # "trust-ncg"
curve       = "wicket3" # "wicket3" # "points2" # "points3"
show_figures = True

if curve == "wicket2":

    r = 1
    l = np.pi*r
    k = 1/r

    print("Setting equality constraints for wicket2 ...")
    print("curve     = ",curve)
    print("radius    = ",r)
    print("length    = ",l)
    print("curvature = ",k)

    bc_dof    = [[0, 1, 3, 4],[ 0, 1, 3,  4]]
    bc_value  = [[r, 0, 0, 1],[-r, 0, 0, -1]]

elif curve == "wicket3":

    r     = 2
    angle = np.pi/2
    k     = 1/r

#    cx    = [r,       0, -r*angle**2, 0, -r*angle,           0]
#    cy    = [0, r*angle,           0, 1,        0, -r*angle**2]


    # 0   1   2   3   4   5
    # x   x'  x'' y   y'  y''
    bc_dof    = [[0, 1, 3, 4],[0,3],[ 0, 1, 3,  4]]
    bc_value  = [[r, 0, 0, 1],[0,r],[-r, 0, 0, -1]]

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
    
    # parameter t in [0,1]
    # x   =       r*cos(k*t)
    # y   =       r*sin(k*t)
    # xt  =    -k*r*sin(k*t)
    # yt  =     k*r*cos(k*t)
    # xtt = -k**2*r*cos(k*t)
    # ytt = -k**2*r*sin(k*t)

    # parameter s in [0,l]
    # x   =    r*cos(s/r)
    # y   =    r*sin(s/r)
    # xt  =     -sin(s/r)
    # yt  =      cos(s/r)
    # xtt = -1/r*cos(s/r)
    # ytt = -1/r*sin(s/r)


    x0[2] = -1/r
    x0[8] =  1/r

elif curve=="wicket3":
    # 0     - 5
    # 6:8   - 9:11
    # 12:14 - 15:17

    # Vertex 1:
    x0[2]  = -1/r

    # Vertex 2:
    x0[7]  = -1
    x0[11] = -1/r

    # Vertex 3:
    x0[14] =  1/r

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
    ds_old = 1
    ds      = [1]
    f = 0
    #----------------------------------------------------#
    # Main loop:
    #----------------------------------------------------#
    f_evals  = []
    ds_evals = []
    for i in range(maxiter):

        #----------------------------------------------------#
        # Compute the gradient:
        #----------------------------------------------------#
        grad_old = grad
        grad     = mvc_integrand_jacobian(x)
        grad     = assign_constraints_grad(grad, bc_dof)

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

        f_old = f
        f     = mvc_objective_function(x)
        f_evals.append(f)

        ds_old = ds[0]
        ds     = arc_length(x)
        ds_evals.append(np.sum(ds))
        #----------------------------------------------------#
        # Check convergence:
        #----------------------------------------------------#
        dds_rel = np.abs(ds[0]-ds_old)/ds_old
        df_rel = np.abs(f-f_old)/f_old
        print(f"Iteration {i:<3} df_rel: {df_rel:<10.4e}  gamma: {gamma:<10.4e}  dds_rel: {dds_rel:<10.4e}")
        if dds_rel < tol:
            break

if show_figures:

    ds = arc_length(x)

    #------------------------------------------------------------#
    # Plotting:
    #------------------------------------------------------------#
    xi = np.arange(0, 1.01, 0.01) 
    H0  = hermite_quintic(0,xi)
    Ht  = hermite_quintic(1,xi)
    Htt = hermite_quintic(2,xi)

    x = np.reshape(x,(len(bc_dof),6))
    #print("x")
    #print(x)


    fig, ax = plt.subplots(2, 2, figsize=(10, 4))
    for i in range(x.shape[0]-1):

        cx = np.zeros(6)
        cx[0:3] = x[i][0:3]
        cx[3:6] = x[i+1][0:3]

        cy = np.zeros(6)
        cy[0:3] = x[i][3:6]
        cy[3:6] = x[i+1][3:6]

#        cx[1] = ds   *cx[1]
#        cx[2] = ds**2*cx[2]
#        cx[4] = ds   *cx[4]
#        cx[5] = ds**2*cx[5]

        cx[1] = ds[i]   *cx[1]
        cx[4] = ds[i]   *cx[4]
        cx[2] = ds[i]**2*cx[2]
        cx[5] = ds[i]**2*cx[5]

        cy[1] = ds[i]   *cy[1]
        cy[4] = ds[i]   *cy[4]
        cy[2] = ds[i]**2*cy[2]
        cy[5] = ds[i]**2*cy[5]

        xx  = H0@cx
        yy  = H0@cy

        cx_t  = Ht@cx
        cy_t  = Ht@cy

        cx_tt  = Htt@cx
        cy_tt  = Htt@cy

        k  = curvature(cx_t, cy_t, cx_tt, cy_tt)
        ss = xi*ds[i] + np.sum(ds[0:i])


        
        ax[0,0].plot(ss, xx)
        ax[0,0].set_xlabel('s')
        ax[0,0].set_ylabel('x')
        ax[0,0].grid(True)

        ax[0,1].plot(ss, yy)
        ax[0,1].set_xlabel('s')
        ax[0,1].set_ylabel('y')
        ax[0,1].grid(True)

        ax[1,0].plot(ss, k)
        ax[1,0].set_xlabel('s')
        ax[1,0].set_ylabel('curvature')
        ax[1,0].grid(True)

#        angle = np.arange(angles[0],angles[2],arc_length/100)
#        xxx,yyy,tmp,tmp,tmp,tmp = circle_arc_param_u(radius,center,angles[0],arc_length,angle)
#        plt.plot(xxx,yyy,'k')


#        for i in range(num_points):
#
#            ax[1,1].plot(points[i,0],points[i,1],'bo')

        ax[1,1].plot(xx, yy)
        ax[1,1].axis('equal')
        ax[1,1].set_xlabel('x')
        ax[1,1].set_ylabel('y')
        ax[1,1].grid(True)

    #plt.figure()
    #plt.plot(f_evals[1:])
#
    #plt.figure()
    #plt.plot(ds_evals[1:])

    plt.show()