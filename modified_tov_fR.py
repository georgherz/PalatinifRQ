import stellarstructureeqs_fR as sse_fR
import sly_ply_fps_new as spf
import numpy as np
import pdb

def modified_TOV(t, yn, para):
    
    #pdb.set_trace()
    
    dm = sse_fR.Mr(para[0], yn[0], t, yn[1], sse_fR.R, sse_fR.Pr0, sse_fR.Pr, sse_fR.alpha, sse_fR.beta, sse_fR.alpha_P, sse_fR.beta_P, sse_fR.tau_rr, sse_fR.tau_tt, sse_fR.tau_rr_P, sse_fR.tau_tt_P, sse_fR.phi, sse_fR.phi_P, sse_fR.f, sse_fR.f_P, sse_fR.fR, sse_fR.fR_r, sse_fR.fR_P, sse_fR.fR_PP, para[1], sse_fR.rho_P, sse_fR.rho_PP, para[3], para[4], sse_fR.A, sse_fR.a, sse_fR.b, sse_fR.d, sse_fR.e, sse_fR.j, sse_fR.k, sse_fR.l, para[2])
    
    dp = sse_fR.Pr(para[0], yn[0], t, yn[1], sse_fR.R, sse_fR.Pr0, sse_fR.alpha, sse_fR.beta, sse_fR.tau_rr, sse_fR.tau_tt, sse_fR.f, sse_fR.fR, sse_fR.fR_P, para[1], sse_fR.rho_P, para[3], para[2])
    
    return np.array([dp,dm])


print("Es ist alles nur alles nur reine Spekulation")
