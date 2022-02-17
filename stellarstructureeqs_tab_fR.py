import numpy as np
import math
import pdb

#TODO: change sign in R, fR_P, f_P, fR_PP

def R(rho, P, kappa):
    
    #pdb.set_trace()
    
    res = kappa*(rho-3.0*P)
    
    return res

def f(rho, P, R, kappa, alpha1): #alpha1 is the parameter which controls the strength of R^2 in f(R)=R+alpha*R^2
    
    res = R(rho, P, kappa)+alpha1*R(rho, P, kappa)*R(rho, P, kappa)
    
    return res

def fR(rho, P, R, kappa, alpha1): #alpha1 is the parameter which controls the strength of R^2 in f(R)=R+alpha*R^2
    
    res = 1.0+2.0*alpha1*R(rho, P, kappa)
    
    return res

def rho_P(rho, P, dEOS):
    
    res = dEOS(P) # dEOS is already my interpolation function which returns my derivative of the density if I feed the pressure into it
    #res = 4000000000000000*res
    
    return res

def rho_PP(rho, P, rho_P, dEOS, d2EOS):
    
    res = d2EOS(P) # d2EOS is already the interpolation function which returns the second derivative of the density if the pressure is fed into it
    #res = 400000000000000000000000000000*res
    
    
    return res

#def rho_P(rho, P, dEOS):    # Derivative of density according to pressure (d/dp rho) as in Sylvester eqs. 126 to 128.
#                            # rho is the density, P the pressure, dEOS is the derivative of the equation of state
#    gsi = np.log10(rho/(7.426*np.power(10.0,-29.0)))
#    
#    rho_int = rho/(7.426*np.power(10.0,-29.0))
#    P_int = P/(8.262*math.pow(10,-50))
#    
#    #pdb.set_trace()
#    
#    #res = rho/(P*dEOS(gsi))
#    res = rho_int/(P_int*dEOS(gsi))
#
#    resi = res*(8.94*math.pow(10,20))
#    
#    #pdb.set_trace()
#    
#    return resi

# def rho_PP(rho, P, rho_P, dEOS, d2EOS):
#    
#    gsi = np.log10(rho/(7.426*np.power(10.0,-29.0)))
#    rho_int = rho/(7.426*np.power(10.0,-29.0))
#    P_int = P/(8.262*math.pow(10,-50))
#    rho_P_int = rho_P(rho, P, dEOS)/(8.94*math.pow(10,20))
#
#    #pdb.set_trace()
#    
#    #res = -(P*dEOS(gsi)*dEOS(gsi)/(rho*rho)-P*dEOS(gsi)/(rho*rho)+P*d2EOS(gsi)/(rho*rho))/(rho_P(rho, P, dEOS)*rho_P(rho, P, dEOS)*rho_P(rho, P, dEOS))
#    res = -(P_int*dEOS(gsi)*dEOS(gsi)/(rho_int*rho_int)-P_int*dEOS(gsi)/(rho_int*rho_int)+P_int*d2EOS(gsi)/(rho_int*rho_int))*(rho_P_int*rho_P_int*rho_P_int)
#    resi = res*(1.07*math.pow(10,70))
#    
#    return resi

def fR_P(rho, P, kappa, rho_P, dEOS, alpha1):
    
    res = 2.0*alpha1*kappa*(rho_P(rho, P, dEOS)-3.0)
    
    return res

def f_P(rho, P, R, fR, kappa, rho_P, dEOS, alpha1):
    
    res = kappa*(rho_P(rho, P, dEOS)-3.0)*fR(rho, P, R, kappa, alpha1)
    
    return res

def fR_PP(rho, P, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1):
    
    res = 2.0*alpha1*kappa*rho_PP(rho, P, rho_P, dEOS, d2EOS)
    
    return res

# TODO: CHECK WHAT HAPPENS IF I ONLY USE ONE KAPPA IN TAU_RR AND TAU_TT AND DERIVATIVES

def tau_tt(rho, P, R, f, fR, kappa, alpha1):
    
    # res = 1.0/np.sqrt(fR(rho, P, R, kappa, alpha1)*fR(rho, P, R, kappa, alpha1)*fR(rho, P, R, kappa, alpha1)*fR(rho, P, R, kappa, alpha1)/16.0)*(f(rho, P, R, kappa, alpha1)/2.0-kappa*rho)
    res = (1.0/(fR(rho, P, R, kappa, alpha1)*fR(rho, P, R, kappa, alpha1)))*(f(rho, P, R, kappa, alpha1)/2.0-kappa*rho)
    
    return res

def tau_rr(rho, P, R, f, fR, kappa, alpha1):
    
    # res = 1.0/np.sqrt(fR(rho, P, R, kappa, alpha1)*fR(rho, P, R, kappa, alpha1)*fR(rho, P, R, kappa, alpha1)*fR(rho, P, R, kappa, alpha1)/16.0)*(f(rho, P, R, kappa, alpha1)/2.0+kappa*P)
    res = (1.0/(fR(rho, P, R, kappa, alpha1)*fR(rho, P, R, kappa, alpha1)))*(f(rho, P, R, kappa, alpha1)/2.0+kappa*P)
    
    return res

def alpha(rho, P, fR, fR_P, kappa, rho_P, dEOS, alpha1):
    
    res = (rho+P)*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)/fR(rho, P, R, kappa, alpha1)
    
    return res

def beta(rho, P, r, fR, fR_P, kappa, rho_P, dEOS, alpha1):
    
    res = 2.0*r*(fR_P(rho, P, kappa, rho_P, dEOS, alpha1)/fR(rho, P, R, kappa, alpha1))*(1.0-(rho+P)*3.0*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)/(4.0*fR(rho, P, R, kappa, alpha1)))
    
    return res

def Pr0(rho, P, r, m, R, tau_rr, tau_tt, f, fR, kappa, alpha1):
    
    res = ((rho+P)/(r*(r-2.0*m)))*(m-(tau_rr(rho, P, R, f, fR, kappa, alpha1)+tau_tt(rho, P, R, f, fR, kappa, alpha1))*r*r*r/4.0)
    
    return res

def Pr(rho, P, r, m, R, Pr0, alpha, beta, tau_rr, tau_tt, f, fR, fR_P, kappa, rho_P, dEOS, alpha1):
    
    res = -2.0*Pr0(rho, P, r, m, R, tau_rr, tau_tt, f, fR, kappa, alpha1)/((1.0-alpha(rho, P, fR, fR_P, kappa, rho_P, dEOS, alpha1))*(1.0+np.sqrt(1.0-beta(rho, P, r, fR, fR_P, kappa, rho_P, dEOS, alpha1)*Pr0(rho, P, r, m, R, tau_rr, tau_tt, f, fR, kappa, alpha1))))
    
    return res

def A(r, m):
    
    res = (1.0-2.0*m/r)
    
    return res

def fR_r(rho, P, r, m, R, Pr0, Pr, alpha, beta, tau_rr, tau_tt, f, fR, fR_P, kappa, rho_P, dEOS, alpha1):
    
    res = fR_P(rho, P, kappa, rho_P, dEOS, alpha1)*Pr(rho, P, r, m, R, Pr0, alpha, beta, tau_rr, tau_tt, f, fR, fR_P, kappa, rho_P, dEOS, alpha1)
    
    return res

def phi(rho, P, R, tau_rr, tau_tt, f, fR, kappa, alpha1):
    
    res = tau_rr(rho, P, R, f, fR, kappa, alpha1)+tau_tt(rho, P, R, f, fR, kappa, alpha1)
    
    return res

def alpha_P(rho, P, R, fR, fR_P, fR_PP, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1):
    
    res = -fR_P(rho, P, kappa, rho_P, dEOS, alpha1)*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)*(P+rho)/(fR(rho, P, R, kappa, alpha1)*fR(rho, P, R, kappa, alpha1))+(P+rho)*fR_PP(rho, P, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1)/fR(rho, P, R, kappa, alpha1)+(fR_P(rho, P, kappa, rho_P, dEOS, alpha1)*(1.0+rho_P(rho, P, dEOS)))/fR(rho, P, R, kappa, alpha1)
    
    return res

def beta_P(rho, P, r, R, fR, fR_P, fR_PP, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1):
    
    #res = -(2.0*r*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)*(1.0-(3.0*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)*(rho+P))/(4.0*fR(rho, P, R, kappa, alpha1))))/(fR(rho, P, R, kappa, alpha1)*fR(rho, P, R, kappa, alpha1))+2.0*r*fR_PP(rho, P, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1)*(1.0-(3.0*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)*(rho+P))/(4.0*fR(rho, P, R, kappa, alpha1)))/fR(rho, P, R, kappa, alpha1)+2.0*r*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)*(3.0*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)*(P+rho)/(4.0*fR(rho, P, R, kappa, alpha1)*fR(rho, P, R, kappa, alpha1))-3.0*fR_PP(rho, P, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1)*(P+rho)/(4.0*fR(rho, P, R, kappa, alpha1))-3.0*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)*(1.0+rho_P(rho, P, dEOS))/(4.0*fR(rho, P, R, kappa, alpha1)))/fR(rho, P, R, kappa, alpha1)
    
    res = r/(2.0*fR(rho, P, R, kappa, alpha1)*fR(rho, P, R, kappa, alpha1)*fR(rho, P, R, kappa, alpha1))*(6.0*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)*(P+rho)+4.0*fR(rho, P, R, kappa, alpha1)*fR(rho, P, R, kappa, alpha1)*fR_PP(rho, P, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1)-fR(rho, P, R, kappa, alpha1)*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)*(4.0*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)+6.0*fR_PP(rho, P, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1)*(P+rho)+3.0*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)*(1.0+rho_P(rho, P, dEOS))))
    
    return res

def tau_rr_P(rho, P, R, f, f_P, fR, fR_P, kappa, rho_P, dEOS, alpha1):
    
    res = (fR(rho, P, R, kappa, alpha1)*(2.0*kappa+f_P(rho, P, R, fR, kappa, rho_P, dEOS, alpha1))-2.0*(2.0*P*kappa+f(rho, P, R, kappa, alpha1))*fR_P(rho, P, kappa, rho_P, dEOS, alpha1))/(2.0*fR(rho, P, R, kappa, alpha1)*fR(rho, P, R, kappa, alpha1)*fR(rho, P, R, kappa, alpha1))
    
    return res

def tau_tt_P(rho, P, R, f, f_P, fR, fR_P, kappa, rho_P, dEOS, alpha1):
    
    res = (-2.0*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)*(f(rho, P, R, kappa, alpha1)-2.0*kappa*rho)+fR(rho, P, R, kappa, alpha1)*(f_P(rho, P, R, fR, kappa, rho_P, dEOS, alpha1)-2.0*kappa*rho_P(rho, P, dEOS)))/(2.0*fR(rho, P, R, kappa, alpha1)*fR(rho, P, R, kappa, alpha1)*fR(rho, P, R, kappa, alpha1))
    
    return res

def phi_P(rho, P, R, tau_rr_P, tau_tt_p, f, f_P, fR, fR_P, kappa, rho_P, dEOS, alpha1):
    
    res = tau_rr_P(rho, P, R, f, f_P, fR, fR_P, kappa, rho_P, dEOS, alpha1)+tau_tt_P(rho, P, R, f, f_P, fR, fR_P, kappa, rho_P, dEOS, alpha1)
    
    return res
    
def a(rho, P, r, m, R, Pr0, beta, tau_rr, tau_tt, f, fR, fR_P, kappa, rho_P, dEOS, alpha1):
    
    res = 1.0+beta(rho, P, r, fR, fR_P, kappa, rho_P, dEOS, alpha1)*Pr0(rho, P, r, m, R, tau_rr, tau_tt, f, fR, kappa, alpha1)/(2.0*np.sqrt(1.0-beta(rho, P, r, fR, fR_P, kappa, rho_P, dEOS, alpha1)*Pr0(rho, P, r, m, R, tau_rr, tau_tt, f, fR, kappa, alpha1))*(1.0+np.sqrt(1.0-beta(rho, P, r, fR, fR_P, kappa, rho_P, dEOS, alpha1)*Pr0(rho, P, r, m, R, tau_rr, tau_tt, f, fR, kappa, alpha1))))
    
    return res

def b(rho, P, r, m, R, Pr, Pr0, alpha, beta, alpha_P, beta_P, tau_rr, tau_tt, f, fR, fR_P, fR_PP, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1):
    
    res = alpha_P(rho, P, R, fR, fR_P, fR_PP, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1)*Pr(rho, P, r, m, R, Pr0, alpha, beta, tau_rr, tau_tt, f, fR, fR_P, kappa, rho_P, dEOS, alpha1)/(1.0-alpha(rho, P, fR, fR_P, kappa, rho_P, dEOS, alpha1))+beta_P(rho, P, r, R, fR, fR_P, fR_PP, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1)*Pr(rho, P, r, m, R, Pr0, alpha, beta, tau_rr, tau_tt, f, fR, fR_P, kappa, rho_P, dEOS, alpha1)*Pr0(rho, P, r, m, R, tau_rr, tau_tt, f, fR, kappa, alpha1)/(2.0*np.sqrt(1.0-beta(rho, P, r, fR, fR_P, kappa, rho_P, dEOS, alpha1)*Pr0(rho, P, r, m, R, tau_rr, tau_tt, f, fR, kappa, alpha1))*(1.0+np.sqrt(1.0-beta(rho, P, r, fR, fR_P, kappa, rho_P, dEOS, alpha1)*Pr0(rho, P, r, m, R, tau_rr, tau_tt, f, fR, kappa, alpha1))))
    
    return res

def e(rho, P, r, m, R, Pr, Pr0, alpha, beta, tau_rr, tau_tt, tau_rr_P, tau_tt_P, phi, phi_P, f, f_P, fR, fR_P, kappa, rho_P, dEOS, alpha1):
    
    res = ((1.0+rho_P(rho, P, dEOS))/(rho+P)-phi_P(rho, P, R, tau_rr_P, tau_tt_P, f, f_P, fR, fR_P, kappa, rho_P, dEOS, alpha1)*r*r*r/(4.0*(m-phi(rho, P, R, tau_rr, tau_tt, f, fR, kappa, alpha1)*r*r*r/4.0)))*Pr(rho, P, r, m, R, Pr0, alpha, beta, tau_rr, tau_tt, f, fR, fR_P, kappa, rho_P, dEOS, alpha1)-(2.0*(r-m)/(r*(r-2.0*m))+3.0*phi(rho, P, R, tau_rr, tau_tt, f, fR, kappa, alpha1)*r*r/(4.0*(m-phi(rho, P, R, tau_rr, tau_tt, f, fR, kappa, alpha1)*r*r*r/4.0)))
    
    return res

def d(rho, P, r, m, R, tau_rr, tau_tt, f, fR, kappa, alpha1):
    
    res = 2.0/(r-2.0*m)+1.0/(m-phi(rho, P, R, tau_rr, tau_tt, f, fR, kappa, alpha1)*r*r*r/4.0)
    
    return res

def l(rho, P, r, m, R, Pr0, Pr, alpha, beta, tau_rr, tau_tt, f, fR, fR_r, fR_P, kappa, rho_P, dEOS, alpha1):
    
    res = (fR_r(rho, P, r, m, R, Pr0, Pr, alpha, beta, tau_rr, tau_tt, f, fR, fR_P, kappa, rho_P, dEOS, alpha1)/fR(rho, P, R, kappa, alpha1)+2.0/r)/r
    
    return res

def j(rho, P, R, tau_rr, tau_tt, f, fR, kappa, alpha1):
    
    res = (3.0*tau_rr(rho, P, R, f, fR, kappa, alpha1)-tau_tt(rho, P, R, f, fR, kappa, alpha1))/2.0
    
    return res

def k(rho, P, r, m, R, Pr0, Pr, alpha, beta, tau_rr, tau_tt, f, fR, fR_r, fR_P, kappa, rho_P, dEOS, alpha1):
    
    res = (fR_r(rho, P, r, m, R, Pr0, Pr, alpha, beta, tau_rr, tau_tt, f, fR, fR_P, kappa, rho_P, dEOS, alpha1)/fR(rho, P, R, kappa, alpha1))*((2.0*r-3.0*m)/(r*(r-2.0*m))-3.0*fR_r(rho, P, r, m, R, Pr0, Pr, alpha, beta, tau_rr, tau_tt, f, fR, fR_P, kappa, rho_P, dEOS, alpha1)/(4.0*fR(rho, P, R, kappa, alpha1)))
    
    return res

def Mr(rho, P, r, m, R, Pr0, Pr, alpha, beta, alpha_P, beta_P, tau_rr, tau_tt, tau_rr_P, tau_tt_P, phi, phi_P, f, f_P, fR, fR_r, fR_P, fR_PP, kappa, rho_P, rho_PP, dEOS, d2EOS, A, a, b, d, e, j, k, l, alpha1):
    
    res = (j(rho, P, R, tau_rr, tau_tt, f, fR, kappa, alpha1)+A(r, m)*((fR_PP(rho, P, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1)*Pr(rho, P, r, m, R, Pr0, alpha, beta, tau_rr, tau_tt, f, fR, fR_P, kappa, rho_P, dEOS, alpha1)*Pr(rho, P, r, m, R, Pr0, alpha, beta, tau_rr, tau_tt, f, fR, fR_P, kappa, rho_P, dEOS, alpha1)+fR_P(rho, P, kappa, rho_P, dEOS, alpha1)*Pr(rho, P, r, m, R, Pr0, alpha, beta, tau_rr, tau_tt, f, fR, fR_P, kappa, rho_P, dEOS, alpha1)*(e(rho, P, r, m, R, Pr, Pr0, alpha, beta, tau_rr, tau_tt, tau_rr_P, tau_tt_P, phi, phi_P, f, f_P, fR, fR_P, kappa, rho_P, dEOS, alpha1)*a(rho, P, r, m, R, Pr0, beta, tau_rr, tau_tt, f, fR, fR_P, kappa, rho_P, dEOS, alpha1)+b(rho, P, r, m, R, Pr, Pr0, alpha, beta, alpha_P, beta_P, tau_rr, tau_tt, f, fR, fR_P, fR_PP, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1)))/fR(rho, P, R, kappa, alpha1)+k(rho, P, r, m, R, Pr0, Pr, alpha, beta, tau_rr, tau_tt, f, fR, fR_r, fR_P, kappa, rho_P, dEOS, alpha1)))/(l(rho, P, r, m, R, Pr0, Pr, alpha, beta, tau_rr, tau_tt, f, fR, fR_r, fR_P, kappa, rho_P, dEOS, alpha1)-A(r, m)*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)*Pr(rho, P, r, m, R, Pr0, alpha, beta, tau_rr, tau_tt, f, fR, fR_P, kappa, rho_P, dEOS, alpha1)*d(rho, P, r, m, R, tau_rr, tau_tt, f, fR, kappa, alpha1)*a(rho, P, r, m, R, Pr0, beta, tau_rr, tau_tt, f, fR, fR_P, kappa, rho_P, dEOS, alpha1)/fR(rho, P, R, kappa, alpha1))

    return res





print("Eh ols")
