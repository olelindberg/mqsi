from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt

from mqsi_util              import assign_constraints,assign_constraints_grad,Constants
from mvc_integrand_jacobian import mvc_integrand_jacobian
from mvc_objective_function import mvc_objective_function
from hermite_quintic        import hermite_quintic
from curvature              import curvature
from curve_metrics          import arc_length
from mqsi_constraints       import mqsi_constraints
from mqsi_initial_conditions import mqsi_initial_conditions

# node dof arragement
# 0 - x
# 1 - x'
# 2 - x''
# 3 - y
# 4 - y'
# 5 - y''

maxiter     = 100
tol         = 1e-6
solver_type = "gradient_descent" # "trust-constr" # "SLSQP" # "L-BFGS-B" # "dogleg" # "trust-ncg"
curve       = "wicket2"
show_figures = True

center          = [5,6]
radius          = 10
angles          = [0.2*np.pi,0.6*np.pi,1.9*np.pi]
bc_dof,bc_value = mqsi_constraints(curve,center,radius,angles)

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
print(f"Setting initial condition for {curve} ...")
x0 = np.zeros(Constants.NODE_DOFS * len(bc_dof))
x0 = assign_constraints(x0, bc_dof, bc_value)
x0 = mqsi_initial_conditions(curve,x0,center,radius,angles)

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
    ls = np.sum(arc_length(x))
    f  = mvc_objective_function(x)
    x_norm_init = np.linalg.norm(x)
    x_norm_old  = x_norm_init

    #----------------------------------------------------#
    # Main loop:
    #----------------------------------------------------#
    f_evals  = []
    ls_evals = []
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

        xx = np.zeros((len(x),3))
        xx[:,0] = x
        x_old = x
        x     = x - gamma*grad

        f_old = f
        f     = mvc_objective_function(x)
        f_evals.append(f)

        ls_old = ls
        ds     = arc_length(x)
        ls = np.sum(ds)
        ls_evals.append(ls)

        #----------------------------------------------------#
        # Check convergence:
        #----------------------------------------------------#
#        print(f"ls_old {ls_old}")
#        print(f"ls     {ls    }")

        ls_rel = np.abs(ls-ls_old)/ls_old
        df_rel = np.abs(f-f_old)/f_old
        print(f"Iteration {i:<3} df_rel: {df_rel:<10.4e}  gamma: {gamma:<10.4e}  ls_rel: {ls_rel:<10.4e} ls: {ls}")

        if ls_rel < tol:
            break

        x_norm     = np.linalg.norm(x)
        dx_norm    = np.abs(x_norm - x_norm_old)/x_norm_init
        x_norm_old = x_norm
#        print(f"Iteration {i:<3}  gamma: {gamma:<10.4e}  dx_norm: {dx_norm:<10.4e}")
#        if dx_norm < tol:
#            print(f"Converged in {i} iterations.")
#            break




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
#        xxx,yyy,tmp,tmp,tmp,tmp = arc_parameter_circle(radius,center,angles[0],arc_length,angle)
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