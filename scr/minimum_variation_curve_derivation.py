from sympy import Function, Symbol
from sympy import *

def ReplaceStrings(f):
    fstr = str(f)
#    fstr = fstr.replace("(u, v)","")
    fstr = fstr.replace("Derivative(","D(")
#    fstr = fstr.replace("1.0*","")
    fstr = fstr.replace("(u, a)","")
    fstr = fstr.replace("D(cx, u)","cx_u")
    fstr = fstr.replace("D(cy, u)","cy_u")
    fstr = fstr.replace("D(cx, a, u)","cx_au")
    fstr = fstr.replace("D(cy, a, u)","cy_au")

#    fstr = fstr.replace("D(S1, v)","S1_v")
#    fstr = fstr.replace("D(S2, u)","S2_u")
#    fstr = fstr.replace("D(S2, v)","S2_v")
#    fstr = fstr.replace("D(S3, u)","S3_u")
#    fstr = fstr.replace("D(S3, v)","S3_v")
#
    fstr = fstr.replace("D(cx, (u, 2))","cx_uu")
    fstr = fstr.replace("D(cy, (u, 2))","cy_uu")

    fstr = fstr.replace("D(cx, (u, 3))","cx_uuu")
    fstr = fstr.replace("D(cy, (u, 3))","cy_uuu")

    fstr = fstr.replace("D(cx, a, (u, 2))","cx_uu")
    fstr = fstr.replace("D(cy, a, (u, 2))","cy_uu")

    fstr = fstr.replace("D(cx, a, (u, 3))","cx_uuu")
    fstr = fstr.replace("D(cy, a, (u, 3))","cy_uuu")


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


u = Symbol('u')

Cx0   = Symbol('Cx0')
Cy0   = Symbol('Cy0')
Cx1   = Symbol('Cx1')
Cy1   = Symbol('Cy1')

Cx0p  = Symbol('Cx0p')
Cy0p  = Symbol('Cy0p')
Cx1p  = Symbol('Cx1p')
Cy1p  = Symbol('Cy1p')

Cx0pp = Symbol('Cx0pp')
Cy0pp = Symbol('Cy0pp')
Cx1pp = Symbol('Cx1pp')
Cy1pp = Symbol('Cy1pp')
a = Symbol('a')


H1  = Function('H1')(u)
H2  = Function('H2')(u)
H3  = Function('H3')(u)
H4  = Function('H4')(u)
H5  = Function('H5')(u)
H6  = Function('H6')(u)

H = Matrix([H1,H2,H3,H4,H5,H6]) 

C0   = Matrix([Cx0,   Cy0  ,0])
C1   = Matrix([Cx1,   Cy1  ,0])
C0p  = Matrix([Cx0p,  Cy0p ,0])
C1p  = Matrix([Cx1p,  Cy1p ,0])
C0pp = Matrix([Cx0pp, Cy0pp,0])
C1pp = Matrix([Cx1pp, Cy1pp,0])

C = Matrix([C0.transpose(),C1.transpose(),C0p.transpose(),C1p.transpose(),C0pp.transpose(),C1pp.transpose()])

c = C.transpose() * H
cx  = Function('cx')(u,a)
cy  = Function('cy')(u,a)
c   = Matrix([cx,cy,0])


cp   = diff(c  , u)
cpp  = diff(cp , u)
cppp = diff(cpp, u)

cpnorm = (cp.dot(cp))**0.5

# curvature
k = cpp.dot(cp) / (cpnorm**3)

u  = cp.cross(cpp)
du = cp.cross(cppp)
v  = cpnorm**3
dv = 3*cpnorm*cp.dot(cpp)

# derivative of curvature
kt = ((v*du - u*dv)/v**2)[2]

# integrand of curvature variation integral
integrand = kt**2 /cpnorm

# partial derivatives
d1  = diff(integrand, a  )

print("integrand = " + ReplaceStrings(integrand) + ";" )
print("  ")
print("d1   = " + ReplaceStrings(d1) + ";" )
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