import stellarstructereqs_tab_fRQ as sse_fRQ
import sly_ply_fps_new as spf
import numpy as np
import pdb

def modified_TOV(t, yn, para):
    
    dm = sse_fRQ.Mr(para[0], yn[0], t, yn[1], sse_fRQ.Lambda, sse_fRQ.Lambda_P, sse_fRQ.Lambda_PP, sse_fRQ.sigma1, sse_fRQ.sigma2, sse_fRQ.sigma1_P, sse_fRQ.sigma2_P, sse_fRQ.sigma1_PP, sse_fRQ.sigma2_PP, sse_fRQ.omega, sse_fRQ.omega_P, sse_fRQ.omega_PP, sse_fRQ.S, sse_fRQ.S_P, sse_fRQ.S_PP, sse_fRQ.Pr0, sse_fRQ.Pr, sse_fRQ.alpha, sse_fRQ.beta, sse_fRQ.alpha_P, sse_fRQ.beta_P, sse_fRQ.tau_rr, sse_fRQ.tau_tt, sse_fRQ.tau_rr_P, sse_fRQ.tau_tt_P, sse_fRQ.phi, sse_fRQ.phi_P, sse_fRQ.f, sse_fRQ.f_P, sse_fRQ.fR, sse_fRQ.f_snake, sse_fRQ.fR_P, sse_fRQ.fR_PP, sse_fRQ.R, sse_fRQ.Q, para[1], sse_fRQ.rho_P, sse_fRQ.rho_PP, para[4], para[5], sse_fRQ.A, sse_fRQ.a, sse_fRQ.b, sse_fRQ.d, sse_fRQ.e, sse_fRQ.j, sse_fRQ.k, sse_fRQ.l, para[2], para[3])
    
    dp = sse_fRQ.Pr(para[0], yn[0], t, yn[1], sse_fRQ.Lambda, sse_fRQ.Lambda_P, sse_fRQ.sigma1, sse_fRQ.sigma2, sse_fRQ.sigma1_P, sse_fRQ.sigma2_P, sse_fRQ.omega, sse_fRQ.omega_P, sse_fRQ.S, sse_fRQ.S_P, sse_fRQ.Pr0, sse_fRQ.alpha, sse_fRQ.beta, sse_fRQ.tau_rr, sse_fRQ.tau_tt, sse_fRQ.f, sse_fRQ.fR, sse_fRQ.f_snake, sse_fRQ.fR_P, sse_fRQ.R, sse_fRQ.Q, para[1], sse_fRQ.rho_P, para[4], para[2], para[3])
    
    return np.array([dp,dm])


print("Es ist alles nur alles nur reine Spekulation")
