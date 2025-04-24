from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt

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

node_dof = 6
maxiter  = 10000
tol     = 1e-4


curve = "wicket2"
if curve == "wicket2":
    r = 1.0
    l = np.pi*r
    k = 1/r
    print("curve     = ",curve)
    print("radius    = ",r)
    print("length    = ",l)
    print("curvature = ",k)



    bc_dof    = [[0, 1, 3,         4],[ 0, 1, 3,          4]]
    bc_value  = [[r, 0, 0, r*np.pi/2],[-r, 0, 0, -r*np.pi/2]]
elif curve == "wicket3":
    bc_dof    = [[   0,   1,   3,   4],[],[  0,   1,   3,    4]]
    bc_value  = [[-1.0, 0.0, 0.0, 1.0],[],[1.0, 0.0, 0.0, -1.0]]
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
        k = i*node_dof + bc_dof[i][j]
        eq = {'type': 'eq', 'fun': lambda x,i=i,j=j,k=k: x[k] - bc_value[i][j]}
        equality_constraints.append(eq)

#-----------------------------------------------------------------------------#
# Initial condition
#-----------------------------------------------------------------------------#
x0 = np.zeros(node_dof * len(bc_dof))
for i in range(len(bc_dof)):
    for j in range(len(bc_dof[i])):
        k = i*node_dof + bc_dof[i][j]
        x0[k] = bc_value[i][j]

if curve=="wicket2":
    
    print("Setting initial condition for wicket2 ...")
    
    x0[2] = -r/4*np.pi**2
    x0[8] =  r/4*np.pi**2

elif curve=="wicket3":
    x0[7] = 1 # derivative at the middle point
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
options = {"maxiter" : maxiter,'disp': True, "verbose" : 1}
method =  'trust-constr' # 'SLSQP' #
res     = minimize(mvc_objective_function, x0, jac=mvc_integrand_jacobian,constraints=equality_constraints , tol=tol, method=method,options=options)

xi = np.arange(-1, 1.01, 0.01)
H0  = hermite_quintic(0,xi)
Ht  = hermite_quintic(1,xi)
Htt = hermite_quintic(2,xi)

x = np.reshape(res.x,(len(bc_dof),6))
print("x")
print(x)

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