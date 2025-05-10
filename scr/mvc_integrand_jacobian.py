import numpy as np
from mvc_integrand_gradient import mvc_integrand_gradient
from mvc_vertex_to_curve    import mvc_vertex_to_curve
from hermite_quintic        import hermite_quintic

def mvc_integrand_jacobian(x,s):

    jac = 0*x

    gauss_xi, gauss_w = np.polynomial.legendre.leggauss(20)


    for i in range(0,len(x)-6,6):

        x = np.reshape(x,(3,6))
        ii = int(i/6)
        cx,cy = mvc_vertex_to_curve(x,ii)
        x = x.flatten()

        s_xi = 0.5*(s[ii+1] - s[ii])

        Ht   = hermite_quintic(1, gauss_xi,s_xi)
        Htt  = hermite_quintic(2, gauss_xi,s_xi)
        Httt = hermite_quintic(3, gauss_xi,s_xi)


        cx_t   = Ht@cx
        cy_t   = Ht@cy
        cx_tt  = Htt@cx
        cy_tt  = Htt@cy
        cx_ttt = Httt@cx
        cy_ttt = Httt@cy

        fx_a = np.zeros(6)
        fy_a = np.zeros(6)
    
        for j in range(0,6):



            cx_at = Ht[:,j]
            cy_at = Ht[:,j]

            cx_att = Htt[:,j]
            cy_att = Htt[:,j]

            cx_attt = Httt[:,j]
            cy_attt = Httt[:,j]

            fx_a[j] = np.sum(mvc_integrand_gradient(cx_t,cy_t,cx_tt,cy_tt,cx_ttt,cy_ttt,cx_at,    0,cx_att,     0,cx_attt,      0)*gauss_w)
            fy_a[j] = np.sum(mvc_integrand_gradient(cx_t,cy_t,cx_tt,cy_tt,cx_ttt,cy_ttt,    0,cy_at,     0,cy_att,      0,cy_attt)*gauss_w)

        jac[i+0:i+3]  = jac[i+0:i+3]  + fx_a[0:3]
        jac[i+6:i+9]  = jac[i+6:i+9]  + fx_a[3:6]
        jac[i+3:i+6]  = jac[i+3:i+6]  + fy_a[0:3]
        jac[i+9:i+12] = jac[i+9:i+12] + fy_a[3:6]


    return jac
