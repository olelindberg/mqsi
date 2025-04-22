from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt

def curvature(dx,dy,ddx,ddy):

    numerator = np.abs(dx * ddy - dy * ddx)
    denominator = (dx**2 + dy**2)**1.5
    
    return numerator / denominator

def hermite_quintic(nderivative, x):
    x = np.asarray(x)
    H = np.zeros((6, len(x)))

    if nderivative == 0:
        H[0, :] = 1/2 - (15*x)/16 + (5*x**3)/8 - (3*x**5)/16
        H[1, :] = 5/16 - (7*x)/16 - (3*x**2)/8 + (5*x**3)/8 + x**4/16 - (3*x**5)/16
        H[2, :] = 1/16 - x/16 - x**2/8 + x**3/8 + x**4/16 - x**5/16
        H[3, :] = 1/2 + (15*x)/16 - (5*x**3)/8 + (3*x**5)/16
        H[4, :] = -5/16 - (7*x)/16 + (3*x**2)/8 + (5*x**3)/8 - x**4/16 - (3*x**5)/16
        H[5, :] = 1/16 + x/16 - x**2/8 - x**3/8 + x**4/16 + x**5/16

    elif nderivative == 1:
        H[0, :] = -15/16 + (15*x**2)/8 - (15*x**4)/16
        H[1, :] = -7/16 - (3*x)/4 + (15*x**2)/8 + x**3/4 - (15*x**4)/16
        H[2, :] = -1/16 - x/4 + (3*x**2)/8 + x**3/4 - (5*x**4)/16
        H[3, :] = 15/16 - (15*x**2)/8 + (15*x**4)/16
        H[4, :] = -7/16 + (3*x)/4 + (15*x**2)/8 - x**3/4 - (15*x**4)/16
        H[5, :] = 1/16 - x/4 - (3*x**2)/8 + x**3/4 + (5*x**4)/16

    elif nderivative == 2:
        H[0, :] = (15*x)/4 - (15*x**3)/4
        H[1, :] = -3/4 + (15*x)/4 + (3*x**2)/4 - (15*x**3)/4
        H[2, :] = -1/4 + (3*x)/4 + (3*x**2)/4 - (5*x**3)/4
        H[3, :] = -15*x/4 + (15*x**3)/4
        H[4, :] = 3/4 + (15*x)/4 - (3*x**2)/4 - (15*x**3)/4
        H[5, :] = -1/4 - (3*x)/4 + (3*x**2)/4 + (5*x**3)/4

    elif nderivative == 3:
        H[0, :] = 15/4 - (45*x**2)/4
        H[1, :] = 15/4 + (6*x)/4 - (45*x**2)/4
        H[2, :] = 3/4 + (6*x)/4 - (15*x**2)/4
        H[3, :] = -15/4 + (45*x**2)/4
        H[4, :] = 15/4 - (6*x)/4 - (45*x**2)/4
        H[5, :] = -3/4 + (6*x)/4 + (15*x**2)/4

    # Transpose to match MATLAB's H = H'
    return H.T



def curvature_derivative(Cx_t,Cy_t,Cx_tt,Cy_tt,Cx_ttt,Cy_ttt):

    cross_term  = Cx_t * Cy_tt - Cy_t * Cx_tt

    numerator_1 = -3.0 * (Cx_t * Cx_tt + Cy_t * Cy_tt)
    denom_1     = (Cx_t**2 + Cy_t**2)**2.5

    numerator_2 = Cx_t * Cy_ttt - Cy_t * Cx_ttt
    denom_2     = (Cx_t**2 + Cy_t**2)**1.5

    Kt = numerator_1 * cross_term / denom_1 + numerator_2 / denom_2

    return Kt


# Define a function of a vector (array)
def objective(x):

    gauss_xi, gauss_w = np.polynomial.legendre.leggauss(20)

    H1 = hermite_quintic(1, gauss_xi)
    H2 = hermite_quintic(2, gauss_xi)
    H3 = hermite_quintic(3, gauss_xi)

    f = 0.0
    for i in range(0,len(x)-6,6):

        cx = np.zeros(6)
        cx[0:3] = x[i+0:i+3]
        cx[3:6] = x[i+6:i+9]

#        print("cx")
#        print(cx)


        cy = np.zeros(6)
        cy[0:3] = x[i+3:i+6]
        cy[3:6] = x[i+9:i+12]
   #     print("cy")
   #     print(cy)

        cx_t   = H1@cx
        cy_t   = H1@cy
        cx_tt  = H2@cx
        cy_tt  = H2@cy
        cx_ttt = H3@cx
        cy_ttt = H3@cy

        k  = curvature(cx_t, cy_t, cx_tt, cy_tt)
        kt = curvature_derivative(cx_t, cy_t, cx_tt, cy_tt, cx_ttt, cy_ttt)

 #       print("kt")
 #       print(kt)
        f = f + np.sum(kt*kt*gauss_w)


    return f

# node dof arragement
# 0 - x
# 1 - x'
# 2 - x''
# 3 - y
# 4 - y'
# 5 - y''



node_dof = 6
nnodes   = 3

ndof = node_dof * nnodes

# Boundary conditions
bc_dof    = [   0,   1,   3,   4,   0,   1,   3,    4]
bc_value  = [-1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, -1.0]

bc_offset_end = int(ndof - node_dof)
bc_offset     = [0,0,0,0,bc_offset_end,bc_offset_end,bc_offset_end,bc_offset_end] 

# Define the equality constraints:
equality_constraints = []
for i in range(len(bc_dof)):
    eq = {'type': 'eq', 'fun': lambda x,i=i: x[bc_dof[i] + bc_offset[i]] - bc_value[i]}
    equality_constraints.append(eq)

#exit()

# Initial guess
x0 = np.zeros(ndof)
for i in range(len(bc_dof)):
    x0[bc_dof[i] + bc_offset[i]] = bc_value[i]



x0[7] = 1
x0[9] = 0.6



tol     = 1e-16
options = {"maxiter" : 1000,'disp': True}
res     = minimize(objective, x0,constraints=equality_constraints , tol=tol, method='SLSQP',options=options)

xi = np.arange(-1, 1.01, 0.01)
H0 = hermite_quintic(0,xi)

x = res.x
#x = x0
for i in range(0,len(x)-6,6):

    cx = np.zeros(6)
    cx[0:3] = x[i+0:i+3]
    cx[3:6] = x[i+6:i+9]

    cy = np.zeros(6)
    cy[0:3] = x[i+3:i+6]
    cy[3:6] = x[i+9:i+12]

    xx  = H0@cx
    yy  = H0@cy

    plt.figure(1)
    plt.plot(xi, xx)
    plt.title('x')
    plt.grid(True)

    plt.figure(2)
    plt.plot(xi, yy)
    plt.title('y')
    plt.grid(True)

    plt.figure(3)
    plt.plot(xx, yy)
    plt.title('x-y')
    plt.grid(True)

plt.show()