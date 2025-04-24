from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt

def curvature_variation_gradient(Cx_t,Cy_t,Cx_tt,Cy_tt,Cx_ttt,Cy_ttt, Cx_ta, Cy_ta,Cx_tta, Cy_tta,Cx_ttta, Cy_ttta):

    norm_Ct_sq = Cx_t**2 + Cy_t**2
    norm_Ct_2p5 = norm_Ct_sq**2.5
    norm_Ct_3p5 = norm_Ct_sq**3.5
    norm_Ct_1p5 = norm_Ct_sq**1.5
    norm_Ct_0p5 = norm_Ct_sq**0.5

    # Curvature derivative Kt (only z-component is non-zero)
    Kz_t = (-3.0 * (Cx_t * Cx_tt + Cy_t * Cy_tt) * (Cx_t * Cy_tt - Cy_t * Cx_tt) / norm_Ct_2p5 +
             (Cx_t * Cy_ttt - Cy_t * Cx_ttt) / norm_Ct_1p5)

    # Curvature derivative Kta (only z-component is non-zero)
    Kz_ta = (15.0 * (Cx_t * Cx_ta + Cy_t * Cy_ta) * (Cx_t * Cx_tt + Cy_t * Cy_tt) * (Cx_t * Cy_tt - Cy_t * Cx_tt) / norm_Ct_3p5
             - 3.0 * (Cx_t * Cx_ta + Cy_t * Cy_ta) * (Cx_t * Cy_ttt - Cy_t * Cx_ttt) / norm_Ct_2p5
             - 3.0 * (Cx_t * Cx_tt + Cy_t * Cy_tt) * (Cx_t * Cy_tta - Cy_t * Cx_tta + Cx_ta * Cy_tt - Cx_tt * Cy_ta) / norm_Ct_2p5
             - 3.0 * (Cx_t * Cy_tt - Cy_t * Cx_tt) * (Cx_t * Cx_tta + Cy_t * Cy_tta + Cx_ta * Cx_tt + Cy_ta * Cy_tt) / norm_Ct_2p5
             + (Cx_t * Cy_ttta - Cy_t * Cx_ttta + Cx_ta * Cy_ttt - Cy_ta * Cx_ttt) / norm_Ct_1p5)

    # Since Kt and Kta only have z-components in MATLAB, we set x and y components as 0
    Kx_t = Ky_t = 0.0
    Kx_ta = Ky_ta = 0.0

    # Compute fa (scalar)
    fa = (- (Cx_t * Cx_ta + Cy_t * Cy_ta) / norm_Ct_1p5 * (Kx_t**2 + Ky_t**2 + Kz_t**2)
          + 1.0 / norm_Ct_0p5 * (2 * Kx_t * Kx_ta + 2 * Ky_t * Ky_ta + 2 * Kz_t * Kz_ta))

    return fa



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

#    print("function")

    gauss_xi, gauss_w = np.polynomial.legendre.leggauss(20)

    Ht = hermite_quintic(1, gauss_xi)
    Htt = hermite_quintic(2, gauss_xi)
    Httt = hermite_quintic(3, gauss_xi)

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

        cx_t   = Ht@cx
        cy_t   = Ht@cy
        cx_tt  = Htt@cx
        cy_tt  = Htt@cy
        cx_ttt = Httt@cx
        cy_ttt = Httt@cy

        kt = curvature_derivative(cx_t, cy_t, cx_tt, cy_tt, cx_ttt, cy_ttt)

 #       print("kt")
 #       print(kt)
        f = f + np.sum(kt*kt*gauss_w)

    return f

def gradient(x):

#    print("gradient")
    grad = 0*x

    gauss_xi, gauss_w = np.polynomial.legendre.leggauss(20)

    Ht   = hermite_quintic(1, gauss_xi)
    Htt  = hermite_quintic(2, gauss_xi)
    Httt = hermite_quintic(3, gauss_xi)

    for i in range(0,len(x)-6,6):

        cx = np.zeros(6)
        cx[0:3] = x[i+0:i+3]
        cx[3:6] = x[i+6:i+9]

        cy = np.zeros(6)
        cy[0:3] = x[i+3:i+6]
        cy[3:6] = x[i+9:i+12]

        cx_t   = Ht@cx
        cy_t   = Ht@cy
        cx_tt  = Htt@cx
        cy_tt  = Htt@cy
        cx_ttt = Httt@cx
        cy_ttt = Httt@cy

        fx_a = np.zeros(6)
        fy_a = np.zeros(6)
        for j in range(0,6):


            cx_ta = Ht[:,j]
            cy_ta = Ht[:,j]

            cx_tta = Htt[:,j]
            cy_tta = Htt[:,j]

            cx_ttta = Httt[:,j]
            cy_ttta = Httt[:,j]

            fx_a[j] = np.sum(curvature_variation_gradient(cx_t,cy_t,cx_tt,cy_tt,cx_ttt,cy_ttt,cx_ta,0,cx_tta,0,cx_ttta,0)*gauss_w)
            fy_a[j] = np.sum(curvature_variation_gradient(cx_t,cy_t,cx_tt,cy_tt,cx_ttt,cy_ttt,0,cy_ta,0,cy_tta,0,cy_ttta)*gauss_w)

        grad[i+0:i+3]  = grad[i+0:i+3]  + fx_a[0:3]
        grad[i+6:i+9]  = grad[i+6:i+9]  + fx_a[3:6]
        grad[i+3:i+6]  = grad[i+3:i+6]  + fy_a[0:3]
        grad[i+9:i+12] = grad[i+9:i+12] + fy_a[3:6]

    return grad


# node dof arragement
# 0 - x
# 1 - x'
# 2 - x''
# 3 - y
# 4 - y'
# 5 - y''



node_dof = 6
maxiter  = 1000
tol     = 1e-8

curve = "wicket3"
curve = "wicket2"
curve = "points2"
curve = "points3"
if curve == "wicket2":
    bc_dof    = [[   0,   1,   3,   4],[  0,   1,   3,    4]]
    bc_value  = [[-1.0, 0.0, 0.0, 1.0],[1.0, 0.0, 0.0, -1.0]]
elif curve == "wicket3":
    bc_dof    = [[   0,   1,   3,   4],[],[  0,   1,   3,    4]]
    bc_value  = [[-1.0, 0.0, 0.0, 1.0],[],[1.0, 0.0, 0.0, -1.0]]
elif curve == "points2":
    bc_dof    = [[  0,   1,   3,   4],[  0,   1,   3,   4]]
    bc_value  = [[0.0, 1.0, 0.0, 0.0],[1.0, 1.0, 1.0, 0.0]]
elif curve == "points3":
    bc_dof    = [[  0,   3],[0,    3],[  0,   3]]
    bc_value  = [[0.0, 0.0],[0.5,0.25],[1.0, 1.0]]


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

if curve=="wicket3":
    x0[7] = 1 # derivative at the middle point
elif curve=="points3":
    x0[7]  = 1 # derivative at the middle point
    x0[10] = 1 # derivative at the middle point

#-----------------------------------------------------------------------------#
# Solve:
#-----------------------------------------------------------------------------#
options = {"maxiter" : maxiter,'disp': True, "verbose" : 2}
method = 'trust-constr' # 'SLSQP'
res     = minimize(objective, x0, jac=gradient,constraints=equality_constraints , tol=tol, method=method,options=options)

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