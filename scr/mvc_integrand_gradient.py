def mvc_integrand_gradient(cx_t,cy_t,cx_tt,cy_tt,cx_ttt,cy_ttt, cx_at, cy_at,cx_att, cy_att,cx_attt, cy_attt):

    v           =  (cx_t**2 + cy_t**2)**1.5
    dv          =  3*(cx_t*cx_tt + cy_t*cy_tt)*(cx_t**2 + cy_t**2)**0.5
    u           =  cx_t*cy_tt - cx_tt*cy_t
    du          =  cx_t*cy_ttt - cx_ttt*cy_t
    cpnorm      =  (cx_t**2 + cy_t**2)**0.5
    cpnorm_a    =  (cx_t*cx_at + cy_t*cy_at)/(cx_t**2 + cy_t**2)**0.5
    v_a         =  (3.0*cx_t*cx_at + 3.0*cy_t*cy_at)*(cx_t**2 + cy_t**2)**0.5
    dv_a        =  3*(cx_t*cx_tt + cy_t*cy_tt)*(cx_t*cx_at + cy_t*cy_at)/(cx_t**2 + cy_t**2)**0.5 + 3*(cx_t**2 + cy_t**2)**0.5*(cx_t*cx_att + cx_tt*cx_at + cy_t*cy_att + cy_tt*cy_at)
    u_a         =  cx_t*cy_att - cx_tt*cy_at - cy_t*cx_att + cy_tt*cx_at
    du_a        =  cx_t*cy_attt - cx_ttt*cy_at - cy_t*cx_attt + cy_ttt*cx_at

    integrand_a =  -4*(du*v - dv*u)**2*v_a/(cpnorm*v**5) - (du*v - dv*u)**2*cpnorm_a/(cpnorm**2*v**4) + (du*v - dv*u)*(2*du*v_a - 2*dv*u_a - 2*u*dv_a + 2*v*du_a)/(cpnorm*v**4)

    return integrand_a


import numpy as np

def df_da(x_ss, x_sss, x_ssa, x_sssa,y_ss, y_sss, y_ssa, y_sssa):
    """
    Computes the derivative w.r.t. 'a' of the squared curvature-like expression:
    ((x_ss*x_sss + y_ss*y_sss) / sqrt(x_ss**2 + y_ss**2))**2

    Parameters:
    - x_ss, x_sss: x''(s,a), x'''(s,a)
    - y_ss, y_sss: y''(s,a), y'''(s,a)
    - x_ssa, x_sssa: ∂x''/∂a, ∂x'''/∂a
    - y_ssa, y_sssa: ∂y''/∂a, ∂y'''/∂a

    Returns:
    - df/da
    """

    u     = x_ss * x_sss + y_ss * y_sss
    du_da = x_ssa * x_sss + x_ss * x_sssa + y_ssa * y_sss + y_ss * y_sssa

    v     = np.sqrt(x_ss**2 + y_ss**2)
    dv_da = (x_ss * x_ssa + y_ss * y_ssa) / v

    df_da_val = 2 * u * (v * du_da - u * dv_da) / v**3
    return df_da_val