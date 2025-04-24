from sympy import Function, Symbol
from sympy import *

def ReplaceStrings(f):
    fstr = str(f)
#    fstr = fstr.replace("(u, v)","")
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


#    fstr = fstr.replace("D(S1, v)","S1_v")
#    fstr = fstr.replace("D(S2, t)","S2_t")
#    fstr = fstr.replace("D(S2, v)","S2_v")
#    fstr = fstr.replace("D(S3, t)","S3_t")
#    fstr = fstr.replace("D(S3, v)","S3_v")
#
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
#
#Cx0   = Symbol('Cx0')
#Cy0   = Symbol('Cy0')
#Cx1   = Symbol('Cx1')
#Cy1   = Symbol('Cy1')
#
#Cx0p  = Symbol('Cx0p')
#Cy0p  = Symbol('Cy0p')
#Cx1p  = Symbol('Cx1p')
#Cy1p  = Symbol('Cy1p')
#
#Cx0pp = Symbol('Cx0pp')
#Cy0pp = Symbol('Cy0pp')
#Cx1pp = Symbol('Cx1pp')
#Cy1pp = Symbol('Cy1pp')
#
#
#
#H1  = Function('H1')(t)
#H2  = Function('H2')(t)
#H3  = Function('H3')(t)
#H4  = Function('H4')(t)
#H5  = Function('H5')(t)
#H6  = Function('H6')(t)
#
#H = Matrix([H1,H2,H3,H4,H5,H6]) 
#
#C0   = Matrix([Cx0,   Cy0  ,0])
#C1   = Matrix([Cx1,   Cy1  ,0])
#C0p  = Matrix([Cx0p,  Cy0p ,0])
#C1p  = Matrix([Cx1p,  Cy1p ,0])
#C0pp = Matrix([Cx0pp, Cy0pp,0])
#C1pp = Matrix([Cx1pp, Cy1pp,0])
#
#C = Matrix([C0.transpose(),C1.transpose(),C0p.transpose(),C1p.transpose(),C0pp.transpose(),C1pp.transpose()])
#
#c = C.transpose() * H
#cx  = Function('cx')(t,a)
#cy  = Function('cy')(t,a)
#c   = Matrix([cx,cy,0])
#
#
#cp   = diff(c  , t)
#cpp  = diff(cp , t)
#cppp = diff(cpp, t)
#
#cpnorm = (cp.dot(cp))**0.5
#
## curvature
#k = cpp.dot(cp) / (cpnorm**3)
#
#u  = cp.cross(cpp)
#du = cp.cross(cppp)
#v  = cpnorm**3
#dv = 3*cpnorm*cp.dot(cpp)

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



# Derivative(cx_t)
# Derivative(cx_tt)
# Derivative(cx_ttt)
# 
# Derivative(cx_at)
# Derivative(cx_att)
# Derivative(cx_attt)
# 
# Derivative(cy_t)
# Derivative(cy_tt)
# Derivative(cy_ttt)
# 
# Derivative(cy_at)
# Derivative(cy_att)
# Derivative(cy_attt)





#print("integrand = " + ReplaceStrings(integrand) + ";" )
#print("  ")
#print("d1   = " + ReplaceStrings(d1) + ";" )
#print("F   = " + ReplaceStrings(F) + ";" )
#print("G   = " + ReplaceStrings(G) + ";" )
#print("E_u = " + ReplaceStrings(Eu) + ";")
#print("E_v = " + ReplaceStrings(Ev) + ";")
#print("F_u = " + ReplaceStrings(Fu) + ";")
#print("F_v = " + ReplaceStrings(Fv) + ";")
#print("G_u = " + ReplaceStrings(Gu) + ";")
#print("G_v = " + ReplaceStrings(Gv) + ";")
#
#
#file = open('MVS_FirstFundamentalFormAndDerivatives.txt', 'w')
#file.write("E   = " + ReplaceStrings(E) + ";"   + '\n')
#file.write("F   = " + ReplaceStrings(F) + ";"   + '\n')
#file.write("G   = " + ReplaceStrings(G) + ";"   + '\n')
#file.write("E_u = " + ReplaceStrings(Eu) + ";"  + '\n')
#file.write("E_v = " + ReplaceStrings(Ev) + ";"  + '\n')
#file.write("F_u = " + ReplaceStrings(Fu) + ";"  + '\n')
#file.write("F_v = " + ReplaceStrings(Fv) + ";"  + '\n')
#file.write("G_u = " + ReplaceStrings(Gu) + ";"  + '\n')
#file.write("G_v = " + ReplaceStrings(Gv) + ";"  + '\n')
#file.close()