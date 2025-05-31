from sympy import Function, Symbol
from sympy import *

def ReplaceStrings(f):
    fstr = str(f)
    fstr = fstr.replace("Derivative(","D(")
    fstr = fstr.replace("1.0*","")
    fstr = fstr.replace("(t, a)","")
    fstr = fstr.replace("D(cx, t)","cx_t")
    fstr = fstr.replace("D(cy, t)","cy_t")
    fstr = fstr.replace("D(cx, a, t)","cx_at")
    fstr = fstr.replace("D(cy, a, t)","cy_at")

    fstr = fstr.replace("D(u, a)","u_a")
    fstr = fstr.replace("D(v, a)","v_a")
    fstr = fstr.replace("D(du, a)","du_a")
    fstr = fstr.replace("D(dv, a)","dv_a")
    fstr = fstr.replace("D(cpnorm, a)","cpnorm_a")


    fstr = fstr.replace("D(cx, (t, 2))","cx_tt")
    fstr = fstr.replace("D(cy, (t, 2))","cy_tt")

    fstr = fstr.replace("D(cx, (t, 3))","cx_ttt")
    fstr = fstr.replace("D(cy, (t, 3))","cy_ttt")

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


t = Symbol('t')
a = Symbol('a')

# Step 1: Define the curvature variation integral
u  = Function('u')(t,a)
v  = Function('v')(t,a)
du = Function('du')(t,a)
dv = Function('dv')(t,a)
cpnorm = Function('cpnorm')(t,a)

# derivative of curvature
kt = ((v*du - u*dv)/v**2)

# integrand of curvature variation integral
integrand = kt**2 /cpnorm

# Step 2: Partial derivatives
integrand_a  = diff(integrand, a  )

cx  = Function('cx')(t,a)
cy  = Function('cy')(t,a)
c   = Matrix([cx,cy,0])

cp     = diff(c  , t)
cpp    = diff(cp , t)
cppp   = diff(cpp, t)

cpnorm = (cp.dot(cp))**0.5
u      = cp.cross(cpp)
du     = cp.cross(cppp)
v      = cpnorm**3
dv     = 3*cpnorm*cp.dot(cpp)

# Derivative(cpnorm(t, a), a)
cpnorm_a = diff(cpnorm, a)

# Derivative(v(t, a), a)
v_a = diff(v, a)

# Derivative(u(t, a), a)
u_a = diff(u, a)

# Derivative(dv(t, a), a)
dv_a = diff(dv, a)

# Derivative(du(t, a), a)
du_a = diff(du, a)

print("v           = ", ReplaceStrings(v))
print("dv          = ", ReplaceStrings(dv))
print("u           = ", ReplaceStrings(u[2]))
print("du          = ", ReplaceStrings(du[2]))
print("cpnorm      = ", ReplaceStrings(cpnorm))
print("integrad    = ", ReplaceStrings(integrand))
print("integrand_a = ", ReplaceStrings(integrand_a))
print("cpnorm_a    = ", ReplaceStrings(cpnorm_a))
print("v_a         = ", ReplaceStrings(v_a))
print("dv_a        = ", ReplaceStrings(dv_a))
print("u_a         = ", ReplaceStrings(u_a[2]))
print("du_a        = ", ReplaceStrings(du_a[2]))