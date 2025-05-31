import numpy as np
from mvc_integrand_gradient import mvc_integrand_gradient
from hermite_quintic        import hermite_quintic
from curve_metrics          import arc_length

def mvc_integrand_jacobian(x):


    print(arc_length(x))
    ds = arc_length(x)[0]

    jac  = 0*x
#    mass = 0*x
    gauss_xi, gauss_w = np.polynomial.legendre.leggauss(20)

    t    = (gauss_xi+1)/2
    Ht   = hermite_quintic(1, t)
    Htt  = hermite_quintic(2, t)
    Httt = hermite_quintic(3, t)

    for i in range(0,len(x)-6,6):

        cx = np.zeros(6)
        cx[0:3] = x[i+0:i+3]
        cx[3:6] = x[i+6:i+9]

        cy = np.zeros(6)
        cy[0:3] = x[i+3:i+6]
        cy[3:6] = x[i+9:i+12]

        cx[1] = ds*cx[1]
        cx[4] = ds*cx[4]
        cx[2] = ds**2*cx[2]
        cx[5] = ds**2*cx[5]


        
        cy[1] = ds*cy[1]
        cy[4] = ds*cy[4]
        cy[2] = ds**2*cy[2]
        cy[5] = ds**2*cy[5]

        cx_t   = Ht@cx
        cy_t   = Ht@cy
        cx_tt  = Htt@cx
        cy_tt  = Htt@cy
        cx_ttt = Httt@cx
        cy_ttt = Httt@cy

        fx_a = np.zeros(6)
        fy_a = np.zeros(6)
        m    = np.zeros(6)
        for j in range(0,6):

            dss = 1
            if j==1 or j==4:
                dss = ds

            if j==2 or j==5:
                dss = ds**2

            cx_at = dss*Ht[:,j]
            cy_at = dss*Ht[:,j]

            cx_att = dss*Htt[:,j]
            cy_att = dss*Htt[:,j]

            cx_attt = dss*Httt[:,j]
            cy_attt = dss*Httt[:,j]

            fx_a[j] = np.sum(mvc_integrand_gradient(cx_t,cy_t,cx_tt,cy_tt,cx_ttt,cy_ttt,cx_at,    0,cx_att,     0,cx_attt,      0)*gauss_w/2)
            fy_a[j] = np.sum(mvc_integrand_gradient(cx_t,cy_t,cx_tt,cy_tt,cx_ttt,cy_ttt,    0,cy_at,     0,cy_att,      0,cy_attt)*gauss_w/2)
            
#            m[j] = np.sum(ds*gauss_w/2)

        jac[i+0:i+3]  = jac[i+0:i+3]  + fx_a[0:3]
        jac[i+6:i+9]  = jac[i+6:i+9]  + fx_a[3:6]
        jac[i+3:i+6]  = jac[i+3:i+6]  + fy_a[0:3]
        jac[i+9:i+12] = jac[i+9:i+12] + fy_a[3:6]

#        mass[i+0:i+3]  = mass[i+0:i+3]  + m[0:3]
#        mass[i+6:i+9]  = mass[i+6:i+9]  + m[3:6]
#        mass[i+3:i+6]  = mass[i+3:i+6]  + m[0:3]
#        mass[i+9:i+12] = mass[i+9:i+12] + m[3:6]

#    jac = jac/mass
    return jac
