def mvc_integrand(cx_t,cy_t,cx_tt,cy_tt,cx_ttt,cy_ttt):

    v           =  (cx_t**2 + cy_t**2)**1.5
    dv          =  3*(cx_t*cx_tt + cy_t*cy_tt)*(cx_t**2 + cy_t**2)**0.5
    u           =  cx_t*cy_tt - cx_tt*cy_t
    du          =  cx_t*cy_ttt - cx_ttt*cy_t
    cpnorm      =  (cx_t**2 + cy_t**2)**0.5
    integrad    =  (du*v - dv*u)**2/(cpnorm*v**4)

    return integrad