import numpy as np
from mvc_integrand_gradient import mvc_integrand_gradient
from mvc_vertex_to_curve    import mvc_vertex_to_curve
from hermite_quintic        import hermite_quintic
from hermite_quintic_derivatives import hermite_quintic_derivatives
from mvc_integrand_gradient import df_da

def mvc_integrand_jacobian(x,s):

    jac = 0*x

    gausss_xi, gauss_w = np.polynomial.legendre.leggauss(20)
    t = 0.5*(gausss_xi+1)

    for i in range(0,len(x)-6,6):

        x = np.reshape(x,(3,6))
        ii = int(i/6)
        cx,cy = mvc_vertex_to_curve(x,ii)
        x = x.flatten()

        s_t = (s[ii+1] - s[ii])        
        xx, yy, x_s, y_s, x_ss, y_ss,x_sss,y_sss = hermite_quintic_derivatives(s_t,t,cx,cy)

        Hss  = s_t**2*hermite_quintic(2, t,s_t)
        Hsss = s_t**3*hermite_quintic(3, t,s_t)

        fx_a = np.zeros(6)
        fy_a = np.zeros(6)
    
        for j in range(0,6):

            x_ssa = Hss[:,j]
            y_ssa = Hss[:,j]

            x_sssa = Hsss[:,j]
            y_sssa = Hsss[:,j]

            fx_a[j] = s_t*0.5*np.sum(df_da(x_ss, x_sss, x_ssa, x_sssa,y_ss, y_sss,     0,      0)*gauss_w)
            fy_a[j] = s_t*0.5*np.sum(df_da(x_ss, x_sss,     0,      0,y_ss, y_sss, y_ssa, y_sssa)*gauss_w)

        jac[i+0:i+3]  = jac[i+0:i+3]  + fx_a[0:3]
        jac[i+6:i+9]  = jac[i+6:i+9]  + fx_a[3:6]
        jac[i+3:i+6]  = jac[i+3:i+6]  + fy_a[0:3]
        jac[i+9:i+12] = jac[i+9:i+12] + fy_a[3:6]


    return jac
