from hermite_quintic import hermite_quintic

def hermite_quintic_derivatives(s_t,t,cx,cy):

    t_s   = 1/s_t

    H0    = hermite_quintic(0,t,s_t)
    H_t   = hermite_quintic(1,t,s_t)
    H_tt  = hermite_quintic(2,t,s_t)
    H_ttt = hermite_quintic(3,t,s_t)

    xx  = H0@cx
    yy  = H0@cy

    cx_s   = t_s*H_t@cx
    cy_s   = t_s*H_t@cy

    cx_ss  = t_s**2*H_tt@cx
    cy_ss  = t_s**2*H_tt@cy

    cx_ss  = t_s**2*H_tt@cx
    cy_ss  = t_s**2*H_tt@cy

    cx_sss = t_s**3*H_ttt@cx
    cy_sss = t_s**3*H_ttt@cy

    return xx, yy, cx_s, cy_s, cx_ss, cy_ss,cx_sss, cy_sss