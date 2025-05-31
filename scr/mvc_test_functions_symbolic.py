from sympy import Function, Symbol
from sympy import *
from sympy import symbols, Function, sin

def ReplaceStrings(f):
    fstr = str(f)
    fstr = fstr.replace("Derivative(" ,"D(")
    fstr = fstr.replace("1.0*"        ,"")
    fstr = fstr.replace("x(s)"        ,"x")
    fstr = fstr.replace("y(s)"        ,"y")
    fstr = fstr.replace("D(x, s)"     ,"x_s")
    fstr = fstr.replace("D(y, s)"     ,"y_s")
    fstr = fstr.replace("D(x, (s, 2))","x_ss")
    fstr = fstr.replace("D(y, (s, 2))","y_ss")
    fstr = fstr.replace("D(x, (s, 3))","x_sss")
    fstr = fstr.replace("D(y, (s, 3))","y_sss")


    fstr = fstr.replace("D(cx, a, t)","cx_at")
    fstr = fstr.replace("D(cy, a, t)","cy_at")

    fstr = fstr.replace("D(u, a)","u_a")
    fstr = fstr.replace("D(v, a)","v_a")
    fstr = fstr.replace("D(du, a)","du_a")
    fstr = fstr.replace("D(dv, a)","dv_a")
    fstr = fstr.replace("D(cpnorm, a)","cpnorm_a")

    fstr = fstr.replace("D(cx, a, (t, 2))","cx_att")
    fstr = fstr.replace("D(cy, a, (t, 2))","cy_att")

    fstr = fstr.replace("D(cx, a, (t, 3))","cx_attt")
    fstr = fstr.replace("D(cy, a, (t, 3))","cy_attt")


#    fstr = fstr.replace("D(S2, (u, 2))","S2_uu")
#    fstr = fstr.replace("D(S3, (u, 2))","S3_uu")
#
#    fstr = fstr.replace("D(S1, u, v)","S1_uv")
#    fstr = fstr.replace("D(S2, u, v)","S2_uv")
#    fstr = fstr.replace("D(S3, u, v)","S3_uv")
#
#    fstr = fstr.replace("D(S1, (v, 2))","S1_vv")
#    fstr = fstr.replace("D(S2, (v, 2))","S2_vv")
#    fstr = fstr.replace("D(S3, (v, 2))","S3_vv")
#
#    fstr = fstr.replace("**",        "^")
#
    return fstr

is_functions = False

a = Symbol('a')
b = Symbol('b')
s = Symbol('s')

y0 = Symbol('y0')
y1 = Symbol('y1')
m0 = Symbol('m0')
m1 = Symbol('m1')
a0 = Symbol('a0')
a1 = Symbol('a1')

test_function_name = "y2"
test_function_name = "sine-cosine"
test_function_name = "hermite-quintic"
test_function_name = "symbolic"

if test_function_name == "symbolic":
    x = Function('x')(s)
    y = Function('y')(s)
elif test_function_name == "y2":
    x = s
    y = s**2
elif test_function_name == "sine-cosine":
    x = sin(s)
    y = cos(s)
elif test_function_name == "hermite-quintic":

    # Normalized parameter t in [0, 1]
    t = (s - a) / (b - a)
    
    # Basis functions
    h00 = 1 - 10*t**3 + 15*t**4 - 6*t**5
    h10 = t - 6*t**3 + 8*t**4 - 3*t**5
    h20 = 0.5*(t**2 - 3*t**3 + 3*t**4 - t**5)
    h01 = 10*t**3 - 15*t**4 + 6*t**5
    h11 = -4*t**3 + 7*t**4 - 3*t**5
    h21 = 0.5*(t**3 - 2*t**4 + t**5)

    # Hermite quintic interpolation
    H = (h00 * y0 +
         h10 * (b - a) * m0 +
         h20 * (b - a)**2 * a0 +
         h01 * y1 +
         h11 * (b - a) * m1 +
         h21 * (b - a)**2 * a1)

    x = s
    y = H

f1 = x

# calcation the arc length
f2 = sqrt(diff(x,s)**2 + diff(y,s)**2)

# calcation the curvature
f3 = sqrt(diff(x,(s,2))**2 + diff(y,(s,2))**2)

# calcation the curvature variation
f4 = diff(f3,s)



print("f1 = ", f1)
print("f2 = ", ReplaceStrings(f2))
print("f3 = ", ReplaceStrings(f3))
print("f4 = ", ReplaceStrings(f4))

F1 = integrate(f1,(s,a,b))
F2 = integrate(f2,(s,a,b))
F3 = integrate(f3,(s,a,b))
F4 = integrate(f4,(s,a,b))

print("F1 = ", ReplaceStrings(F1))
print("F2 = ", ReplaceStrings(F2))
print("F3 = ", ReplaceStrings(F3))
print("F4 = ", ReplaceStrings(F4))

#plot_parametric((x, y), (s, a, b))