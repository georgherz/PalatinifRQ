import numpy as np
import sly_ply_fps_new as spf
import math
import pdb

#################################################################
#                                                               #
# Stellar Structure Equations for f(R,Q) = R+alpha*R^2+beta*Q   #
# Disclaimer: In all equations I replaced f_Q by beta           #
#             since df(R,Q)/dQ=beta                             #
#                                                               #
#################################################################

#TODO: ACHTUNG: R_p und R_pp hab ich nicht als Funktionen programmiert sondern so reingeschrieben, dh falls ich ein anderes Modell mit einem anderen R nehme muss ich auch R_p und R_pp aendern!
# TODO: check all equations whether they are right, check all the derivatives in mathematica whether they are right (I might have made a mistake back then when I was coding that), find right Q
# TODO: Check all functions whether they have the right dependences later on when using them in other functions (I am changing quite a lot right now)
# TODO: Check in other paper (bouncing cosmologies) whether all the signs are right

def R(rho, P, kappa):
    
    res = kappa*(rho-3.0*P)
    
    return res

def fR(rho, P, R, kappa, alpha1): #alpha1 is the parameter which controls the strength of R^2 in f(R)=R+alpha*R^2+beta*Q. Beta1 controls the strength of Q
    
    res = 1.0+2.0*alpha1*R(rho, P, kappa)
    
    return res

def f_snake(rho, P, fR, R, kappa, alpha1): # f_snake = R+alpha*R^2, i.e. it is the f from f(R)
    
    res = R(rho, P, kappa)+alpha1*R(rho, P, kappa)*R(rho, P, kappa)
    
    return res

def Q(rho, P, fR, f_snake, R, kappa, alpha1, beta1): # Q as defined by eq. 2.16 in Olmo et al. 2012 (unreadable paper). Use - sign where +/- in equation (because it is given that way in Barragan and Olmo 2010 (Isotropic and anisotropic bouncing cosmologies in Palatini gravity)). IMPORTANT: This is actually beta*Q, i.e. I do not need a beta in front of Q in the equation for f.
    
    # Q as in Olmo
    res = -(f_snake(rho, P, fR, R, kappa, alpha1)+fR(rho, P, R, kappa, alpha1)*fR(rho, P, R, kappa, alpha1)/(4.0*beta1)+2.0*kappa*P)+(beta1/16.0)*((3.0*(R(rho, P, kappa)+fR(rho, P, R, kappa, alpha1)/beta1)-np.sqrt((R(rho, P, kappa)+fR(rho, P, R, kappa, alpha1)/beta1)*(R(rho, P, kappa)+fR(rho, P, R, kappa, alpha1)/beta1)-4.0*kappa*(rho+P)/beta1))*(3.0*(R(rho, P, kappa)+fR(rho, P, R, kappa, alpha1)/beta1)-np.sqrt((R(rho, P, kappa)+fR(rho, P, R, kappa, alpha1)/beta1)*(R(rho, P, kappa)+fR(rho, P, R, kappa, alpha1)/beta1)-4.0*kappa*(rho+P)/beta1)))
    
    #pdb.set_trace()
    
    #res = -(f_snake(rho, P, fR, R, kappa, alpha1)+fR(rho, P, R, kappa, alpha1)*fR(rho, P, R, kappa, alpha1)/(4.0*beta1)+2.0*kappa*P)+(1.0/(16.0*beta1))*((3.0*(R(rho, P, kappa)*beta1+fR(rho, P, R, kappa, alpha1))-np.sqrt((R(rho, P, kappa)*beta1+fR(rho, P, R, kappa, alpha1))*(R(rho, P, kappa)*beta1+fR(rho, P, R, kappa, alpha1))-4.0*kappa*(rho+P)*beta1))*(3.0*(R(rho, P, kappa)*beta1+fR(rho, P, R, kappa, alpha1))-np.sqrt((R(rho, P, kappa)*beta1+fR(rho, P, R, kappa, alpha1))*(R(rho, P, kappa)*beta1+fR(rho, P, R, kappa, alpha1))-4.0*kappa*(rho+P)*beta1)))
    
    return res

def f(rho, P, fR, f_snake, R, Q, kappa, alpha1, beta1): #alpha1 is the parameter which controls the strength of R^2 in f(R)=R+alpha*R^2+beta*Q. beta1 controls the strengt of Q
    
    res = R(rho, P, kappa)+alpha1*R(rho, P, kappa)*R(rho, P, kappa)+Q(rho, P, fR, f_snake, R, kappa, alpha1, beta1) # IMPORTANT: Q as defined above is actually beta*Q, i.e. I do not need a beta in front of Q here.
    
    return res

def rho_P(rho, P, dEOS):
    
    res = dEOS(P) # dEOS is already my interpolation function which returns my derivative of the density if I feed the pressure into it
    
    return res

def rho_PP(rho, P, rho_P, dEOS, d2EOS):
    
    res = d2EOS(P) # d2EOS is already the interpolation function which returns the second derivative of the density if the pressure is fed into it
    
    return res

#def rho_P(rho, P, dEOS):    # Derivative of density according to pressure (d/dp rho) as in Sylvester eqs. 126 to 128.
#                            # rho is the density, P the pressure, dEOS is the derivative of the equation of state
#    gsi = np.log10(rho/(7.426*np.power(10.0,-29.0)))
#    rho_int = rho/(7.426*np.power(10.0,-29.0))
#    P_int = P/(8.262*math.pow(10,-50))
#    res = rho_int/(P_int*dEOS(gsi))
#    resi = res*(8.987*math.pow(10,20))
#    
#    return resi
#
#def rho_PP(rho, P, rho_P, dEOS, d2EOS):
#    
#    gsi = np.log10(rho/(7.426*np.power(10.0,-29.0)))
#    rho_int = rho/(7.426*np.power(10.0,-29.0))
#    P_int = P/(8.262*math.pow(10,-50))
#    rho_P_int = rho_P(rho, P, dEOS)/(8.987*math.pow(10,20))
#    #res = -(P*dEOS(gsi)*dEOS(gsi)/(rho*rho)-P*dEOS(gsi)/(rho*rho)+P*d2EOS(gsi)/(rho*rho))/(rho_P(rho, P, dEOS)*rho_P(rho, P, dEOS)*rho_P(rho, P, dEOS))
#    #res = -(P_int*dEOS(gsi)*dEOS(gsi)/(rho_int*rho_int)-P_int*dEOS(gsi)/(rho_int*rho_int)+P_int*d2EOS(gsi)/(rho_int*rho_int))*(rho_P_int*rho_P_int*rho_P_int)
#    res = -(P_int*dEOS(gsi)*dEOS(gsi)/(rho_int*rho_int)-P_int*dEOS(gsi)/(rho_int*rho_int)+P_int*d2EOS(gsi)/(rho_int*rho_int))/(P_int*P_int*P_int*dEOS(gsi)*dEOS(gsi)*dEOS(gsi)/#(rho_int*rho_int*rho_int))
#    resi = res*(1.087*math.pow(10,70))
#    
#    return resi

def fR_P(rho, P, kappa, rho_P, dEOS, alpha1):
    
    res = 2.0*alpha1*kappa*(rho_P(rho, P, dEOS)-3.0)
    
    return res

def fR_PP(rho, P, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1):
    
    res = 2.0*alpha1*kappa*rho_PP(rho, P, rho_P, dEOS, d2EOS)
    
    return res

def f_P(rho, P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1):
    
    #res = kappa*(rho_P(rho, P, dEOS)-3.0)*fR(rho, P, R, kappa, alpha1)+((1.0/8.0)*beta1*(3.0*(fR(rho, P, R, kappa, alpha1)/beta1+R(rho, P, kappa))-np.sqrt(((fR(rho, P, R, kappa, alpha1)+beta1*R(rho, P, kappa))*(fR(rho, P, R, kappa, alpha1)+beta1*R(rho, P, kappa))-4.0*beta1*kappa*(P+rho))/(beta1*beta1)))*(3.0*(fR_P(rho, P, kappa, rho_P, dEOS, alpha1)/beta1+kappa*(rho_P(rho, P, dEOS)-3.0))+(2.0*beta1*kappa*(1.0+rho_P(rho, P, dEOS))-(fR_P(rho, P, kappa, rho_P, dEOS, alpha1)+beta1*kappa*(rho_P(rho, P, dEOS)-3.0))*(fR(rho, P, R, kappa, alpha1)+beta1*R(rho, P, kappa)))/(beta1*beta1*np.sqrt(((fR(rho, P, R, kappa, alpha1)+beta1*R(rho, P, kappa))*(fR(rho, P, R, kappa, alpha1)+beta1*R(rho, P, kappa))-4.0*beta1*kappa*(P+rho))/(beta1*beta1))))-2.0*kappa-fR(rho, P, R, kappa, alpha1)*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)/(2.0*beta1)-kappa*(rho_P(rho, P, dEOS)-3.0)*fR(rho, P, R, kappa, alpha1))
    
    # Without simplify
    res = kappa*(rho_P(rho, P, dEOS)-3.0)*fR(rho, P, R, kappa, alpha1)+((beta1/8.0)*(3.0*(fR(rho, P, R, kappa, alpha1)/beta1+R(rho, P, kappa))-np.sqrt((fR(rho, P, R, kappa, alpha1)/beta1+R(rho, P, kappa))*(fR(rho, P, R, kappa, alpha1)/beta1+R(rho, P, kappa))-4.0*kappa*(P+rho)/beta1))*(3.0*(fR_P(rho, P, kappa, rho_P, dEOS, alpha1)/beta1+kappa*(rho_P(rho, P, dEOS)-3.0))-(2.0*(fR_P(rho, P, kappa, rho_P, dEOS, alpha1)/beta1+kappa*(rho_P(rho, P, dEOS)-3.0))*(fR(rho, P, R, kappa, alpha1)/beta1+R(rho, P, kappa))-4.0*kappa*(1.0+rho_P(rho, P, dEOS))/beta1)/(2.0*np.sqrt((fR(rho, P, R, kappa, alpha1)/beta1+R(rho, P, kappa))*(fR(rho, P, R, kappa, alpha1)/beta1+R(rho, P, kappa))-4.0*kappa*(P+rho)/beta1)))-2.0*kappa-fR(rho, P, R, kappa, alpha1)*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)/(2.0*beta1)-kappa*(rho_P(rho, P, dEOS)-3.0)*fR(rho, P, R, kappa, alpha1))
    
    #res = kappa*(rho_P(rho, P, dEOS)-3.0)*fR(rho, P, R, kappa, alpha1)+((1.0/(beta1*8.0))*(3.0*(fR(rho, P, R, kappa, alpha1)+beta1*R(rho, P, kappa))-np.sqrt((fR(rho, P, R, kappa, alpha1)+beta1*R(rho, P, kappa))*(fR(rho, P, R, kappa, alpha1)+beta1*R(rho, P, kappa))-4.0*kappa*(P+rho)*beta1))*(3.0*beta1*kappa*(rho_P(rho, P, dEOS)-3.0)-(2.0*beta1*kappa*(rho_P(rho, P, dEOS)-3.0)*(fR(rho, P, R, kappa, alpha1)+beta1*R(rho, P, kappa))-4.0*kappa*(1.0+rho_P(rho, P, dEOS))*beta1)/(2.0*np.sqrt((fR(rho, P, R, kappa, alpha1)+beta1*R(rho, P, kappa))*(fR(rho, P, R, kappa, alpha1)+beta1*R(rho, P, kappa))-4.0*kappa*(P+rho)*beta1)))-2.0*kappa-fR(rho, P, R, kappa, alpha1)*fR_P(rho, P, kappa, rho_P, dEOS, alpha1)/(2.0*beta1)-kappa*(rho_P(rho, P, dEOS)-3.0)*fR(rho, P, R, kappa, alpha1))
    
    #pdb.set_trace()
    #print(str(res1)+"    "+str(res2)+"    "+str(res))
    
    return res

# TODO: Ab hier die Gleichungen weiter korrigieren

def Lambda(rho, P, fR, R, kappa, alpha1, beta1): # definition of my function Lambda as defined in Olmo et al. 2012 (eq. 2.15) --> fQ=beta
    
    res = (np.sqrt(2.0*beta1)/8.0)*(3.0*(R(rho, P, kappa)+fR(rho, P, R, kappa, alpha1)/beta1)-np.sqrt((R(rho, P, kappa)+fR(rho, P, R, kappa, alpha1)/beta1)*(R(rho, P, kappa)+fR(rho, P, R, kappa, alpha1)/beta1)-4.0*kappa*(rho+P)/beta1))
    
    return res

def sigma1(rho, P, Lambda, fR, R, kappa, alpha1, beta1): # definition of sigma1 (eq. 2.12 Olmo e al.)
    
    res = fR(rho, P, R, kappa, alpha1)/2.0+np.sqrt(2.0*beta1)*np.sqrt(Lambda(rho, P, fR, R, kappa, alpha1, beta1)*Lambda(rho, P, fR, R, kappa, alpha1, beta1)-kappa*(rho+P))
    
    return res

def sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1): # definition of sigma1 (eq. 2.12 Olmo et al. 2012)
    
    res = fR(rho, P, R, kappa, alpha1)/2.0+np.sqrt(2.0*beta1)*Lambda(rho, P, fR, R, kappa, alpha1, beta1)
    
    return res

def omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1): # definition of omega (eq. 3.2 Olmo et al. 2012)
    
    res = np.sqrt(sigma1(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1))
    
    return res

def S(rho, P, Lambda, sigma1, sigma2, omega, fR, R, kappa, alpha1, beta1): # definition of S (eq 3.1 Olmo et al. 2012)
    
    res = sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)/omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)
    
    return res

def Lambda_P(rho, P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1): # derivative of lambda (eq. 2.15 Olmo et al. 2012). Derivative was calculate using mathematica
    
    res = np.sqrt(beta1/32.0)*(3.0*(fR_P(rho, P, kappa, rho_P, dEOS, alpha1)/beta1+kappa*(rho_P(rho, P, dEOS)-3.0))-(2.0*(fR(rho, P, R, kappa, alpha1)/beta1+R(rho, P, kappa))*(fR_P(rho, P, kappa, rho_P, dEOS, alpha1)/beta1+kappa*(rho_P(rho, P, dEOS)-3.0))-4.0*kappa*(1.0+rho_P(rho, P, dEOS))/beta1)/(2.0*np.sqrt((fR(rho, P, R, kappa, alpha1)/beta1+R(rho, P, kappa))*(fR(rho, P, R, kappa, alpha1)/beta1+R(rho, P, kappa))-4.0*kappa*(P+rho)/beta1)))
    
    return res

def sigma1_P(rho, P, Lambda, Lambda_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1):
    
    res = fR_P(rho, P, kappa, rho_P, dEOS, alpha1)/2.0+np.sqrt(beta1/2.0)*(2.0*Lambda(rho, P, fR, R, kappa, alpha1, beta1)*Lambda_P(rho, P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)-kappa*(1.0+rho_P(rho, P, dEOS)))/np.sqrt(Lambda(rho, P, fR, R, kappa, alpha1, beta1)*Lambda(rho, P, fR, R, kappa, alpha1, beta1)-kappa*(P+rho))
    
    return res

def sigma2_P(rho, P, Lambda_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1):
    
    res = fR_P(rho, P, kappa, rho_P, dEOS, alpha1)/2.0+np.sqrt(2.0*beta1)*Lambda_P(rho, P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)
    
    return res

def omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1):
    
    res = 0.5*(sigma1_P(rho, P, Lambda, Lambda_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)+sigma1(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2_P(rho, P, Lambda_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1))/omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)
    
    return res


def S_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1): # first derivative of 3.1 (olmo et al. 2012), derivative has been calculated with mathematica
    
    res = (2.0*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2_P(rho, P, Lambda_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)-sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1))/(omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)*omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1))
    
    return res

def Lambda_PP(rho, P, fR, fR_P, fR_PP, R, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1): # Second derivative of Lambda. Derivative was taken with mathematica
    
    res = np.sqrt(beta1/32.0)*(3.0*(fR_PP(rho, P, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1)/beta1+kappa*rho_PP(rho, P, rho_P, dEOS, d2EOS))+((2.0*(fR(rho, P, R, kappa, alpha1)/beta1+R(rho, P, kappa))*(fR_P(rho, P, kappa, rho_P, dEOS, alpha1)/beta1+kappa*(rho_P(rho, P, dEOS)-3.0))-4.0*kappa*(1.0+rho_P(rho, P, dEOS))/beta1)*(2.0*(fR(rho, P, R, kappa, alpha1)/beta1+R(rho, P, kappa))*(fR_P(rho, P, kappa, rho_P, dEOS, alpha1)/beta1+kappa*(rho_P(rho, P, dEOS)-3.0))-4.0*kappa*(1.0+rho_P(rho, P, dEOS))/beta1))/(4.0*((fR(rho, P, R, kappa, alpha1)/beta1+R(rho, P, kappa))*(fR(rho, P, R, kappa, alpha1)/beta1+R(rho, P, kappa))-4.0*kappa*(P+rho)/beta1)*np.sqrt((fR(rho, P, R, kappa, alpha1)/beta1+R(rho, P, kappa))*(fR(rho, P, R, kappa, alpha1)/beta1+R(rho, P, kappa))-4.0*kappa*(P+rho)/beta1))-(2.0*(fR_P(rho, P, kappa, rho_P, dEOS, alpha1)/beta1+kappa*(rho_P(rho, P, dEOS)-3.0))*(fR_P(rho, P, kappa, rho_P, dEOS, alpha1)/beta1+kappa*(rho_P(rho, P, dEOS)-3.0))+2.0*(fR(rho, P, R, kappa, alpha1)/beta1+R(rho, P, kappa))*(fR_PP(rho, P, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1)/beta1+kappa*rho_PP(rho, P, rho_P, dEOS, d2EOS))-4.0*kappa*rho_PP(rho, P, rho_P, dEOS, d2EOS)/beta1)/(2.0*np.sqrt((fR(rho, P, R, kappa, alpha1)/beta1+R(rho, P, kappa))*(fR(rho, P, R, kappa, alpha1)/beta1+R(rho, P, kappa))-4.0*kappa*(P+rho)/beta1)))
        
    return res

def sigma1_PP(rho, P, Lambda, Lambda_P, Lambda_PP, fR, fR_P, fR_PP, R, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1):
    
    res = fR_PP(rho, P, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1)/2.0-(np.sqrt(2.0*beta1)*(2.0*Lambda(rho, P, fR, R, kappa, alpha1, beta1)*Lambda_P(rho, P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)-kappa*(1.0+rho_P(rho, P, dEOS)))*(2.0*Lambda(rho, P, fR, R, kappa, alpha1, beta1)*Lambda_P(rho, P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)-kappa*(1.0+rho_P(rho, P, dEOS))))/(4.0*(Lambda(rho, P, fR, R, kappa, alpha1, beta1)*Lambda(rho, P, fR, R, kappa, alpha1, beta1)-kappa*(P+rho))*np.sqrt((Lambda(rho, P, fR, R, kappa, alpha1, beta1)*Lambda(rho, P, fR, R, kappa, alpha1, beta1)-kappa*(P+rho))))+(np.sqrt(2.0*beta1)*(2.0*Lambda_P(rho, P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*Lambda_P(rho, P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)+2.0*Lambda(rho, P, fR, R, kappa, alpha1, beta1)*Lambda_PP(rho, P, fR, fR_P, fR_PP, R, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1)-kappa*rho_PP(rho, P, rho_P, dEOS, d2EOS)))/(2.0*np.sqrt((Lambda(rho, P, fR, R, kappa, alpha1, beta1)*Lambda(rho, P, fR, R, kappa, alpha1, beta1)-kappa*(P+rho))))
    
    return res

def sigma2_PP(rho, P, Lambda_PP, fR, fR_P, fR_PP, R, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1):
    
    res = fR_PP(rho, P, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1)/2.0+np.sqrt(2.0*beta1)*Lambda_PP(rho, P, fR, fR_P, fR_PP, R, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1)
    
    return res

def omega_PP(rho, P, Lambda, Lambda_P, Lambda_PP, sigma1, sigma2, sigma1_P, sigma2_P, sigma1_PP, sigma2_PP, omega, omega_P, fR, fR_P, fR_PP, R, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1):
    
    res = -omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)+(sigma1_PP(rho, P, Lambda, Lambda_P, Lambda_PP, fR, fR_P, fR_PP, R, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)+2.0*sigma1_P(rho, P, Lambda, Lambda_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*sigma2_P(rho, P, Lambda_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)+sigma1(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2_PP(rho, P, Lambda_PP, fR, fR_P, fR_PP, R, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1))/(2.0*omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1))

    return res

def S_PP(rho, P, Lambda, Lambda_P, Lambda_PP, sigma1, sigma2, sigma1_P, sigma2_P, sigma1_PP, sigma2_PP, omega, omega_P, omega_PP, fR, fR_P, fR_PP, R, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1): # second derivative of eq. 3.1 from unreadable paper (Olmo et al. 2012), second derivative has been calculated with mathematica
    
    res = (2.0*omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)*omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)*sigma2_P(rho, P, Lambda_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*sigma2_P(rho, P, Lambda_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)+sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*(2.0*omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)-omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)*omega_PP(rho, P, Lambda, Lambda_P, Lambda_PP, sigma1, sigma2, sigma1_P, sigma2_P, sigma1_PP, sigma2_PP, omega, omega_P, fR, fR_P, fR_PP, R, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1))+2.0*omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*(omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)*sigma2_PP(rho, P, Lambda_PP, fR, fR_P, fR_PP, R, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1)-2.0*omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*sigma2_P(rho, P, Lambda_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)))/(omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)*omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)*omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1))
    
    return res

def alpha(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1): # Eq. 3.27 unreadable Paper (Olmo et al. 2012)
    
    res = 0.5*(rho+P)*(omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)+S_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/S(rho, P, Lambda, sigma1, sigma2, omega, fR, R, kappa, alpha1, beta1))
    
    return res

def beta(rho, P, r, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1): # Eq. 3.28 unreadable Paper (Olmo et al. 2012)
    
    res = 2.0*r*(omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1))*(1.0-0.5*(rho+P)*(1.5*omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)-(omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)-S_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/S(rho, P, Lambda, sigma1, sigma2, omega, fR, R, kappa, alpha1, beta1))))
    
    #print(res)
    
    return res

def tau_tt(rho, P, Lambda, sigma1, sigma2, f, fR, f_snake, R, Q, kappa, alpha1, beta1): # Definition from eq. 2.5 (Olmo et al. 2012) and from comment on p. 7 first paragraph in Olmo et al. 2012
    
    res = 1.0/np.sqrt(sigma1(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1))*(f(rho, P, fR, f_snake, R, Q, kappa, alpha1, beta1)/2.0-kappa*rho)
    
    return res

def tau_rr(rho, P, Lambda, sigma1, sigma2, f, fR, f_snake, R, Q, kappa, alpha1, beta1): # Definition from eq. 2.5 (Olmo et al. 2012) and from comment on p. 7 first paragraph in Olmo et al. 2012
    
    res = 1.0/np.sqrt(sigma1(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1))*(f(rho, P, fR, f_snake, R, Q, kappa, alpha1, beta1)/2.0+kappa*P)
    
    return res

def Pr0(rho, P, r, m, Lambda, sigma1, sigma2, omega, S, tau_rr, tau_tt, f, fR, f_snake, R, Q, kappa, alpha1, beta1): # eq. 3.26 from Olmo et al. 2012
    
    res = (rho+P)/(r*(r-2.0*m))*(m-(tau_rr(rho, P, Lambda, sigma1, sigma2, f, fR, f_snake, R, Q, kappa, alpha1, beta1)+omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)*tau_tt(rho, P, Lambda, sigma1, sigma2, f, fR, f_snake, R, Q, kappa, alpha1, beta1)/S(rho, P, Lambda, sigma1, sigma2, omega, fR, R, kappa, alpha1, beta1))*r*r*r/4.0)
    
    return res

def Pr(rho, P, r, m, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, Pr0, alpha, beta, tau_rr, tau_tt, f, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1): # eq. 3.25 from Olmo et al. 2012
    
    #pdb.set_trace()
    
    res = -2.0*Pr0(rho, P, r, m, Lambda, sigma1, sigma2, omega, S, tau_rr, tau_tt, f, fR, f_snake, R, Q, kappa, alpha1, beta1)/((1.0-alpha(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1))*(1.0+np.sqrt(1.0-beta(rho, P, r, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*Pr0(rho, P, r, m, Lambda, sigma1, sigma2, omega, S, tau_rr, tau_tt, f, fR, f_snake, R, Q, kappa, alpha1, beta1))))
    
    return res

def alpha_P(rho, P, Lambda, Lambda_P, Lambda_PP, sigma1, sigma2, sigma1_P, sigma2_P, sigma1_PP, sigma2_PP, omega, omega_P, omega_PP, S, S_P, S_PP, fR, fR_P, fR_PP, R, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1):
    
    res = 0.5*(S_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/S(rho, P, Lambda, sigma1, sigma2, omega, fR, R, kappa, alpha1, beta1)+omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1))*(1.0+rho_P(rho, P, dEOS))+0.5*(P+rho)*(-S_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*S_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/(S(rho, P, Lambda, sigma1, sigma2, omega, fR, R, kappa, alpha1, beta1)*S(rho, P, Lambda, sigma1, sigma2, omega, fR, R, kappa, alpha1, beta1))+S_PP(rho, P, Lambda, Lambda_P, Lambda_PP, sigma1, sigma2, sigma1_P, sigma2_P, sigma1_PP, sigma2_PP, omega, omega_P, omega_PP, fR, fR_P, fR_PP, R, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1)/S(rho, P, Lambda, sigma1, sigma2, omega, fR, R, kappa, alpha1, beta1)-omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/(omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)*omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1))+omega_PP(rho, P, Lambda, Lambda_P, Lambda_PP, sigma1, sigma2, sigma1_P, sigma2_P, sigma1_PP, sigma2_PP, omega, omega_P, fR, fR_P, fR_PP, R, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1)/omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1))
    
    return res

def beta_P(rho, P, r, Lambda, Lambda_P, Lambda_PP, sigma1, sigma2, sigma1_P, sigma2_P, sigma1_PP, sigma2_PP, omega, omega_P, omega_PP, S, S_P, S_PP, fR, fR_P, fR_PP, R, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1):
    
    res = -(2.0*r*omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*(1.0-0.5*(P+rho)*(S_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/S(rho, P, Lambda, sigma1, sigma2, omega, fR, R, kappa, alpha1, beta1)+0.5*omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)))*omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1))/(omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)*omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1))+2.0*r*(1.0-0.5*(P+rho)*(S_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/S(rho, P, Lambda, sigma1, sigma2, omega, fR, R, kappa, alpha1, beta1)+0.5*omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)))*omega_PP(rho, P, Lambda, Lambda_P, Lambda_PP, sigma1, sigma2, sigma1_P, sigma2_P, sigma1_PP, sigma2_PP, omega, omega_P, fR, fR_P, fR_PP, R, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1)/omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)+(2.0*r*omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1))*(-0.5*(S_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/S(rho, P, Lambda, sigma1, sigma2, omega, fR, R, kappa, alpha1, beta1)+0.5*omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1))*(1.0+rho_P(rho, P, dEOS))-0.5*(P+rho)*(-S_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*S_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/(S(rho, P, Lambda, sigma1, sigma2, omega, fR, R, kappa, alpha1, beta1)*S(rho, P, Lambda, sigma1, sigma2, omega, fR, R, kappa, alpha1, beta1))+S_PP(rho, P, Lambda, Lambda_P, Lambda_PP, sigma1, sigma2, sigma1_P, sigma2_P, sigma1_PP, sigma2_PP, omega, omega_P, omega_PP, fR, fR_P, fR_PP, R, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1)/S(rho, P, Lambda, sigma1, sigma2, omega, fR, R, kappa, alpha1, beta1)-omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/(2.0*omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)*omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1))+0.5*omega_PP(rho, P, Lambda, Lambda_P, Lambda_PP, sigma1, sigma2, sigma1_P, sigma2_P, sigma1_PP, sigma2_PP, omega, omega_P, fR, fR_P, fR_PP, R, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1)/omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)))
    
    return res

def tau_tt_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, f, f_P, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1):
    
    res = (f_P(rho, P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/2.0-kappa*rho_P(rho, P, dEOS))/np.sqrt(sigma1(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1))-(f(rho, P, fR, f_snake, R, Q, kappa, alpha1, beta1)/2.0-kappa*rho)*(sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma1_P(rho, P, Lambda, Lambda_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)+3.0*sigma1(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2_P(rho, P, Lambda_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1))/(2.0*sigma1(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*np.sqrt(sigma1(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)))
    
    return res

def tau_rr_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, f, f_P, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1):
    
    res = (kappa+f_P(rho, P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/2.0)/np.sqrt(sigma1(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1))-(f(rho, P, fR, f_snake, R, Q, kappa, alpha1, beta1)/2.0+P*kappa)*(sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma1_P(rho, P, Lambda, Lambda_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)+3.0*sigma1(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2_P(rho, P, Lambda_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1))/(2.0*sigma1(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*np.sqrt(sigma1(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)*sigma2(rho, P, Lambda, fR, R, kappa, alpha1, beta1)))
    
    return res

def phi(rho, P, Lambda, sigma1, sigma2, omega, S, tau_rr, tau_tt, f, fR, f_snake, R, Q, kappa, alpha1, beta1): # as defined on p. 9 of Olmo et al. 2012
    
    res = tau_rr(rho, P, Lambda, sigma1, sigma2, f, fR, f_snake, R, Q, kappa, alpha1, beta1)+omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)*tau_tt(rho, P, Lambda, sigma1, sigma2, f, fR, f_snake, R, Q, kappa, alpha1, beta1)/S(rho, P, Lambda, sigma1, sigma2, omega, fR, R, kappa, alpha1, beta1)
    
    return res

def phi_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, tau_tt, tau_rr_P, tau_tt_P, f, f_P, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1): # derivadive of phi taken with mathematica
    
    res = -tau_tt(rho, P, Lambda, sigma1, sigma2, f, fR, f_snake, R, Q, kappa, alpha1, beta1)*omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)*S_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/(S(rho, P, Lambda, sigma1, sigma2, omega, fR, R, kappa, alpha1, beta1)*S(rho, P, Lambda, sigma1, sigma2, omega, fR, R, kappa, alpha1, beta1))+tau_rr_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, f, f_P, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1)+omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)*tau_tt_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, f, f_P, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1)/S(rho, P, Lambda, sigma1, sigma2, omega, fR, R, kappa, alpha1, beta1)+tau_tt(rho, P, Lambda, sigma1, sigma2, f, fR, f_snake, R, Q, kappa, alpha1, beta1)*omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)/S(rho, P, Lambda, sigma1, sigma2, omega, fR, R, kappa, alpha1, beta1)
    
    return res

def A(r, m): # A(r) is the A one finds in the Schwarzschild Metric (also given on p. 10 of Olmo et al. 2012)
    
    res = (1.0-2.0*m/r)
    
    return res

def a(rho, P, r, m, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, Pr0, beta, tau_rr, tau_tt, f, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1): # definition of a() from Sylvester 2018
    
    res = 1.0+beta(rho, P, r, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*Pr0(rho, P, r, m, Lambda, sigma1, sigma2, omega, S, tau_rr, tau_tt, f, fR, f_snake, R, Q, kappa, alpha1, beta1)/(2.0*np.sqrt(1.0-beta(rho, P, r, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*Pr0(rho, P, r, m, Lambda, sigma1, sigma2, omega, S, tau_rr, tau_tt, f, fR, f_snake, R, Q, kappa, alpha1, beta1))*(1.0+np.sqrt(1.0-beta(rho, P, r, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*Pr0(rho, P, r, m, Lambda, sigma1, sigma2, omega, S, tau_rr, tau_tt, f, fR, f_snake, R, Q, kappa, alpha1, beta1))))
    
    return res

def b(rho, P, r, m, Lambda, Lambda_P, Lambda_PP, sigma1, sigma2, sigma1_P, sigma2_P, sigma1_PP, sigma2_PP, omega, omega_P, omega_PP, S, S_P, S_PP, Pr0, Pr, alpha, beta, alpha_P, beta_P, tau_rr, tau_tt, f, fR, f_snake, fR_P, fR_PP, R, Q, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1): # definition of b() from Sylvester 2018
    
    res = alpha_P(rho, P, Lambda, Lambda_P, Lambda_PP, sigma1, sigma2, sigma1_P, sigma2_P, sigma1_PP, sigma2_PP, omega, omega_P, omega_PP, S, S_P, S_PP, fR, fR_P, fR_PP, R, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1)*Pr(rho, P, r, m, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, Pr0, alpha, beta, tau_rr, tau_tt, f, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1)/(1.0-alpha(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1))+beta_P(rho, P, r, Lambda, Lambda_P, Lambda_PP, sigma1, sigma2, sigma1_P, sigma2_P, sigma1_PP, sigma2_PP, omega, omega_P, omega_PP, S, S_P, S_PP, fR, fR_P, fR_PP, R, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1)*Pr(rho, P, r, m, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, Pr0, alpha, beta, tau_rr, tau_tt, f, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1)*Pr0(rho, P, r, m, Lambda, sigma1, sigma2, omega, S, tau_rr, tau_tt, f, fR, f_snake, R, Q, kappa, alpha1, beta1)/(2.0*np.sqrt(1.0-beta(rho, P, r, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*Pr0(rho, P, r, m, Lambda, sigma1, sigma2, omega, S, tau_rr, tau_tt, f, fR, f_snake, R, Q, kappa, alpha1, beta1))*(1.0+np.sqrt(1.0-beta(rho, P, r, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*Pr0(rho, P, r, m, Lambda, sigma1, sigma2, omega, S, tau_rr, tau_tt, f, fR, f_snake, R, Q, kappa, alpha1, beta1))))
    
    return res

def e(rho, P, r, m, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, Pr0, Pr, alpha, beta, tau_rr, tau_tt, tau_rr_P, tau_tt_P, phi, phi_P, f, f_P, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1):
    
    res = ((1.0+rho_P(rho, P, dEOS))/(rho+P)-phi_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, tau_tt, tau_rr_P, tau_tt_P, f, f_P, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1)*np.power(r,3)/(4.0*(m-phi(rho, P, Lambda, sigma1, sigma2, omega, S, tau_rr, tau_tt, f, fR, f_snake, R, Q, kappa, alpha1, beta1)*np.power(r,3)/4.0)))*Pr(rho, P, r, m, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, Pr0, alpha, beta, tau_rr, tau_tt, f, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1)-(2.0*(r-m)/(r*(r-2.0*m))+3.0*phi(rho, P, Lambda, sigma1, sigma2, omega, S, tau_rr, tau_tt, f, fR, f_snake, R, Q, kappa, alpha1, beta1)*np.power(r,2)/(4.0*(m-phi(rho, P, Lambda, sigma1, sigma2, omega, S, tau_rr, tau_tt, f, fR, f_snake, R, Q, kappa, alpha1, beta1)*np.power(r,3)/4.0)))
    
    return res

def d(rho, P, r, m, Lambda, sigma1, sigma2, omega, S, tau_rr, tau_tt, phi, f, fR, f_snake, R, Q, kappa, alpha1, beta1):
    
    res = 2.0/(r-2.0*m)+1.0/(m-phi(rho, P, Lambda, sigma1, sigma2, omega, S, tau_rr, tau_tt, f, fR, f_snake, R, Q, kappa, alpha1, beta1)*np.power(r,3)/4.0)
    
    return res

def l(rho, P, r, m, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, Pr0, Pr, alpha, beta, tau_rr, tau_tt, f, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1):
    
    res = (omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*Pr(rho, P, r, m, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, Pr0, alpha, beta, tau_rr, tau_tt, f, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1)/omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)+2.0/r)/r
    
    return res

def j(rho, P, Lambda, sigma1, sigma2, omega, S, tau_rr, tau_tt, f, fR, f_snake, R, Q, kappa, alpha1, beta1):
    
    res = (3.0*tau_rr(rho, P, Lambda, sigma1, sigma2, f, fR, f_snake, R, Q, kappa, alpha1, beta1)-omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)*tau_tt(rho, P, Lambda, sigma1, sigma2, f, fR, f_snake, R, Q, kappa, alpha1, beta1)/S(rho, P, Lambda, sigma1, sigma2, omega, fR, R, kappa, alpha1, beta1))/2.0
    
    return res

def k(rho, P, r, m, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, Pr0, Pr, alpha, beta, tau_rr, tau_tt, f, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1):
    
    res = omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*Pr(rho, P, r, m, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, Pr0, alpha, beta, tau_rr, tau_tt, f, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1)*((2.0*r-3.0*m)/(r*(r-2.0*m))-3.0*omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*Pr(rho, P, r, m, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, Pr0, alpha, beta, tau_rr, tau_tt, f, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1)/(4.0*omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)))/omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)
    
    return res

def Mr(rho, P, r, m, Lambda, Lambda_P, Lambda_PP, sigma1, sigma2, sigma1_P, sigma2_P, sigma1_PP, sigma2_PP, omega, omega_P, omega_PP, S, S_P, S_PP, Pr0, Pr, alpha, beta, alpha_P, beta_P, tau_rr, tau_tt, tau_rr_P, tau_tt_P, phi, phi_P, f, f_P, fR, f_snake, fR_P, fR_PP, R, Q, kappa, rho_P, rho_PP, dEOS, d2EOS, A, a, b, d, e, j, k, l, alpha1, beta1):
    
    res = (j(rho, P, Lambda, sigma1, sigma2, omega, S, tau_rr, tau_tt, f, fR, f_snake, R, Q, kappa, alpha1, beta1)+A(r, m)*((omega_PP(rho, P, Lambda, Lambda_P, Lambda_PP, sigma1, sigma2, sigma1_P, sigma2_P, sigma1_PP, sigma2_PP, omega, omega_P, fR, fR_P, fR_PP, R, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1)*Pr(rho, P, r, m, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, Pr0, alpha, beta, tau_rr, tau_tt, f, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1)*Pr(rho, P, r, m, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, Pr0, alpha, beta, tau_rr, tau_tt, f, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1)+omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*Pr(rho, P, r, m, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, Pr0, alpha, beta, tau_rr, tau_tt, f, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1)*(e(rho, P, r, m, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, Pr0, Pr, alpha, beta, tau_rr, tau_tt, tau_rr_P, tau_tt_P, phi, phi_P, f, f_P, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1)*a(rho, P, r, m, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, Pr0, beta, tau_rr, tau_tt, f, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1)+b(rho, P, r, m, Lambda, Lambda_P, Lambda_PP, sigma1, sigma2, sigma1_P, sigma2_P, sigma1_PP, sigma2_PP, omega, omega_P, omega_PP, S, S_P, S_PP, Pr0, Pr, alpha, beta, alpha_P, beta_P, tau_rr, tau_tt, f, fR, f_snake, fR_P, fR_PP, R, Q, kappa, rho_P, rho_PP, dEOS, d2EOS, alpha1, beta1)))/omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1)+k(rho, P, r, m, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, Pr0, Pr, alpha, beta, tau_rr, tau_tt, f, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1)))/(l(rho, P, r, m, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, Pr0, Pr, alpha, beta, tau_rr, tau_tt, f, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1)-A(r, m)*omega_P(rho, P, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, fR, fR_P, R, kappa, rho_P, dEOS, alpha1, beta1)*Pr(rho, P, r, m, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, Pr0, alpha, beta, tau_rr, tau_tt, f, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1)*d(rho, P, r, m, Lambda, sigma1, sigma2, omega, S, tau_rr, tau_tt, phi, f, fR, f_snake, R, Q, kappa, alpha1, beta1)*a(rho, P, r, m, Lambda, Lambda_P, sigma1, sigma2, sigma1_P, sigma2_P, omega, omega_P, S, S_P, Pr0, beta, tau_rr, tau_tt, f, fR, f_snake, fR_P, R, Q, kappa, rho_P, dEOS, alpha1, beta1)/omega(rho, P, Lambda, sigma1, sigma2, fR, R, kappa, alpha1, beta1))
    
    return res



# TODO: check all functions whether I typed them right
# TODO: Might have to change A (and everything that depends on it).
