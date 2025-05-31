import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

test_function_name = "sine-cosine"
test_function_name = "y2"

xi, w = np.polynomial.legendre.leggauss(24)

t = (xi+1)/2
w = w/2

s    = t
s_t  = 1

if test_function_name == "y2":

    x     = s
    y     = s**2
    x_s   = 1
    y_s   = 2*s
    x_ss  = 0
    y_ss  = 2
    x_sss = 0
    y_sss = 0

    F1_exact =  -1/3 + 2*np.sqrt(2)/3
    F2_exact =  np.arcsinh(2)/4 + np.sqrt(5)/2
    F3_exact =  2
    F4_exact =  0

elif test_function_name == "sine-cosine":

    x     = np.sin(s)
    y     = np.cos(s)
    x_s   = np.cos(s)
    y_s   = -np.sin(s)
    x_ss  = -np.sin(s)
    y_ss  = -np.cos(s)
    x_sss = -np.cos(s)
    y_sss = np.sin(s)

    F1_exact =  1.00000000000000
    F2_exact =  1.00000000000000
    F3_exact =  1.00000000000000
    F4_exact =  0

f1 =  np.sqrt(x**2 + y**2)
f2 =  np.sqrt(x_s**2 + y_s**2)
f3 =  np.sqrt(x_ss**2 + y_ss**2)
f4 =  (x_ss*x_sss + y_ss*y_sss)/np.sqrt(x_ss**2 + y_ss**2)

F1 = np.sum(f1*w)*s_t
F2 = np.sum(f2*w)*s_t
F3 = np.sum(f3*w)*s_t
F4 = np.sum(f4*w)*s_t


print("F1_exact   = ", F1_exact)
print("F1_numeric = ", F1)

print("F2_exact   = ", F2_exact)
print("F2_numeric = ", F2)

print("F3_exact   = ", F3_exact)
print("F3_numeric = ", F3)

print("F4_exact   = ", F4_exact)
print("F4_numeric = ", F4)

plt.figure(1)
plt.plot(x, y)
plt.show()

