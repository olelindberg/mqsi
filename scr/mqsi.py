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

    H0 = hermite_quintic(0, gauss_xi)
    H1 = hermite_quintic(1, gauss_xi)
    H2 = hermite_quintic(2, gauss_xi)
    H3 = hermite_quintic(3, gauss_xi)

    #for i in range(0,len(x),12):
    i = 0
    cx = x[i:i+6]
    cy = x[i+6:i+12]

    cx_t   = H1@cx
    cy_t   = H1@cy
    cx_tt  = H2@cx
    cy_tt  = H2@cy
    cx_ttt = H3@cx
    cy_ttt = H3@cy

    k  = curvature(cx_t, cy_t, cx_tt, cy_tt)
    kt = curvature_derivative(cx_t, cy_t, cx_tt, cy_tt, cx_ttt, cy_ttt)


    f = np.sum(kt*kt*gauss_w)

    print(f)

    return f


eq1 = {'type': 'eq', 'fun': lambda x: x[0]  - (0.0)} # x(0)   =  -1
eq2 = {'type': 'eq', 'fun': lambda x: x[1]  - 0.0} # x_t(0) =  0
eq3 = {'type': 'eq', 'fun': lambda x: x[3]  - 1.0} # x(1)   =  1
eq4 = {'type': 'eq', 'fun': lambda x: x[4]  - 0.0} # x_t(1) =  0
eq5 = {'type': 'eq', 'fun': lambda x: x[6]  - 0.0} # y(0)   =  0
eq6 = {'type': 'eq', 'fun': lambda x: x[7]  - 1.0} # y_t(0) =  1
eq7 = {'type': 'eq', 'fun': lambda x: x[9]  - 1.0} # y(1)   =  0
eq8 = {'type': 'eq', 'fun': lambda x: x[10] - (1.0)} # y_t(1) = -1


ndof = 12

x0 = np.zeros(ndof)

x0[0]  = -1.0 # x(0)   = -1
x0[1]  = 0.0 # x_t(0) = 0
x0[3]  = 1.0 # x(1)   = 1
x0[4]  = 0.0 # x_t(1) = 0
x0[6]  = 0.0 # y(0)   = 0
x0[7]  = 1.0 # y_t(0) = 1
x0[9]  = 0.0 # y(1)   = 0
x0[10] = 1.0 # y_t(1) = -1


tol     = 1e-16
options = {"maxiter" : 100,'disp': True}
res     = minimize(objective, x0,constraints=[eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8], tol=tol, method='SLSQP',options=options)

cx = res.x[0:6]
cy = res.x[6:12]

xi = np.arange(-1, 1.01, 0.01)
H0 = hermite_quintic(0,xi)

x  = H0@cx
y  = H0@cy


plt.figure()
plt.plot(xi, x, 'r-')
plt.title('x')

plt.figure()
plt.plot(xi, y, 'r-')
plt.title('y')

plt.figure()
plt.plot(x, y, 'r-')
plt.title('x-y')

plt.grid(True)
plt.show()