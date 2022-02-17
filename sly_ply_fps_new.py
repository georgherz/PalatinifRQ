import numpy as np
import math
import matplotlib.pylab as plt


# definition of the equations of my toy model PLY as it can be found in Pannia et al. (2017)
# Structure of compact stars in R-squared palatini gravity
# returns the value of zeta, dzeta/dgsi and d2zeta/dgsi2 for a given gsi as input
# gsi = log10(rho) --> is the logarithm of the density
# zeta = log10(p) --> is the logarithm of the pressure
def ply_zeta(gsi):
    zeta = 2.0*gsi+5.29355  # calculates value for zeta for given gsi
    return zeta

def ply_dzeta(gsi):
    dzeta = 2.0             # value for dzeta for given gsi is constant 2
    return dzeta

def ply_d2zeta(gsi):
    d2zeta = 0              # value for d2zeta for given gsi is constant zero
    return d2zeta



# defining the equations for the model sly
# the equations for zeta, dzeta and d2zeta are defined separately
def sly_zeta(gsi):
    # defining all the coefficients for the equations
    # as they can be found in Haensel & Potekhin (2004)
    # Analytical representation of unified equations of
    # state of neutron-star matter
    a1 = 6.22
    a2 = 6.121
    a3 = 0.005925
    a4 = 0.16326
    a5 = 6.48
    a6 = 11.4971
    a7 = 19.105
    a8 = 0.8938
    a9 = 6.54
    a10 = 11.4950
    a11 = -22.775
    a12 = 1.5707
    a13 = 4.3
    a14 = 14.08
    a15 = 27.80
    a16 = -1.653
    a17 = 1.50
    a18 = 14.67
    
    # defining the function f0
    # returns the value of f0 at x
    def f0(x):
        
        f = 1.0/(math.exp(x)+1)
        
        return f
    
    zeta = (a1+a2*gsi+a3*math.pow(gsi,3))/(1.0+a4*gsi)*f0(a5*(gsi-a6))+(a7+a8*gsi)*f0(a9*(a10-gsi))+(a11+a12*gsi)*f0(a13*(a14-gsi))+(a15+a16*gsi)*f0(a17*(a18-gsi))
    
    return zeta

def sly_dzeta(gsi):
    a1 = 6.22
    a2 = 6.121
    a3 = 0.005925
    a4 = 0.16326
    a5 = 6.48
    a6 = 11.4971
    a7 = 19.105
    a8 = 0.8938
    a9 = 6.54
    a10 = 11.4950
    a11 = -22.775
    a12 = 1.5707
    a13 = 4.3
    a14 = 14.08
    a15 = 27.80
    a16 = -1.653
    a17 = 1.50
    a18 = 14.67
    
    # defining the function f0
    # returns the value of f0 at x
    def f0(x):
        
        f = 1.0/(math.exp(x)+1)
        
        return f
    
    dzeta = a8*f0(a9*(a10-gsi))+a12*f0(a13*(a14-gsi))+a16*f0(a17*(a18-gsi))+(a13*math.exp(a13*(a14-gsi))*(a11+a12*gsi))*math.pow(f0(a13*(a14-gsi)),2)+(a17*math.exp(a17*(a18-gsi))*(a15+a16*gsi))*math.pow(f0(a17*(a18-gsi)),2)+(a9*math.exp(a9*(a10-gsi))*(a7+a8*gsi))*math.pow(f0(a9*(a10-gsi)),2)+(a2+3.0*a3*math.pow(gsi,2))*f0(a5*(gsi-a6))/(1+a4*gsi)-(a4*(a1+a2*gsi+a3*math.pow(gsi,3)))*f0(a5*(gsi-a6))/math.pow((1.0+a4*gsi),2)-(a5*math.exp(a5*(gsi-a6))*(a1+a2*gsi+a3*math.pow(gsi,3)))*math.pow(f0(a5*(gsi-a6)),2)/(1.0+a4*gsi)
    
    return dzeta

def sly_d2zeta(gsi):
    a1 = 6.22
    a2 = 6.121
    a3 = 0.005925
    a4 = 0.16326
    a5 = 6.48
    a6 = 11.4971
    a7 = 19.105
    a8 = 0.8938
    a9 = 6.54
    a10 = 11.4950
    a11 = -22.775
    a12 = 1.5707
    a13 = 4.3
    a14 = 14.08
    a15 = 27.80
    a16 = -1.653
    a17 = 1.50
    a18 = 14.67
    
    # defining the function f0
    # returns the value of f0 at x
    def f0(x):
        
        f = 1.0/(math.exp(x)+1)
        
        return f
    
    ddzeta = (2.0*a8*a9*math.exp(a9*(a10-gsi)))*math.pow(f0(a9*(a10-gsi)),2)+(2.0*a12*a13*math.exp(a13*(a14-gsi)))*math.pow(f0(a13*(a14-gsi)),2)+(2.0*a16*a17*math.exp(a17*(a18-gsi)))*math.pow(f0(a17*(a18-gsi)),2)+(2.0*math.pow(a13,2)*math.exp(2*a13*(a14-gsi))*(a11+a12*gsi))*math.pow(f0(a13*(a14-gsi)),3)-(math.pow(a13,2)*math.exp(a13*(a14-gsi))*(a11+a12*gsi))*math.pow(f0(a13*(a14-gsi)),2)+(2.0*math.pow(a17,2)*math.exp(2.0*a17*(a18-gsi))*(a15+a16*gsi))*math.pow(f0(a17*(a18-gsi)),3)-(math.pow(a17,2)*math.exp(a17*(a18-gsi))*(a15+a16*gsi))*math.pow(f0(a17*(a18-gsi)),2)+(6.0*a3*gsi)*f0(a5*(-a6+gsi))/(1.0+a4*gsi)+(2.0*math.pow(a9,2)*math.exp(2.0*a9*(a10-gsi))*(a7+a8*gsi))*math.pow(f0(a9*(a10-gsi)),3)-(math.pow(a9,2)*math.exp(a9*(a10-gsi))*(a7+a8*gsi))*math.pow(f0(a9*(a10-gsi)),2)-(2.0*a4*(a2+3.0*a3*math.pow(gsi,2)))*f0(a5*(-a6+gsi))/math.pow((1.0+a4*gsi),2)-(2.0*a5*math.exp(a5*(-a6+gsi))*(a2+3.0*a3*math.pow(gsi,2)))*math.pow(f0(a5*(-a6+gsi)),2)/(1+a4*gsi)+(2.0*math.pow(a4,2)*(a1+a2*gsi+a3*math.pow(gsi,3)))*f0(a5*(-a6+gsi))/math.pow((1.0+a4*gsi),3)+(2.0*a4*a5*math.exp(a5*(-a6+gsi))*(a1+a2*gsi+a3*math.pow(gsi,3)))*math.pow(f0(a5*(-a6+gsi)),2)/math.pow((1.0+a4*gsi),2)+(2.0*math.pow(a5,2)*math.exp(2.0*a5*(-a6+gsi))*(a1+a2*gsi+a3*math.pow(gsi,3)))*math.pow(f0(a5*(-a6+gsi)),3)/(1.0+a4*gsi)-(math.pow(a5,2)*math.exp(a5*(-a6+gsi))*(a1+a2*gsi+a3*math.pow(gsi,3)))*math.pow(f0(a5*(-a6+gsi)),2)/(1.0+a4*gsi)
    
    return ddzeta
    
    
    


# defines the equations for the FPS model
# returns values 
def fps_zeta(gsi):
    
    a1 = 6.22
    a2 = 6.121
    a3 = 0.006004
    a4 = 0.16345
    a5 = 6.50
    a6 = 11.8440
    a7 = 17.24
    a8 = 1.065
    a9 = 6.54
    a10 = 11.8421
    a11 = -22.003
    a12 = 1.5552
    a13 = 9.3
    a14 = 14.19
    a15 = 23.73
    a16 = -1.508
    a17 = 1.79
    a18 = 15.13
    
    # defining the function f0
    # returns the value of f0 at x
    def f0(x):
        
        f = 1.0/(math.exp(x)+1)
        
        return f
    
    zeta = (a1+a2*gsi+a3*math.pow(gsi,3))/(1.0+a4*gsi)*f0(a5*(gsi-a6))+(a7+a8*gsi)*f0(a9*(a10-gsi))+(a11+a12*gsi)*f0(a13*(a14-gsi))+(a15+a16*gsi)*f0(a17*(a18-gsi))
    
    return zeta

def fps_dzeta(gsi):
    
    a1 = 6.22
    a2 = 6.121
    a3 = 0.006004
    a4 = 0.16345
    a5 = 6.50
    a6 = 11.8440
    a7 = 17.24
    a8 = 1.065
    a9 = 6.54
    a10 = 11.8421
    a11 = -22.003
    a12 = 1.5552
    a13 = 9.3
    a14 = 14.19
    a15 = 23.73
    a16 = -1.508
    a17 = 1.79
    a18 = 15.13
    
    # defining the function f0
    # returns the value of f0 at x
    def f0(x):
        
        f = 1.0/(math.exp(x)+1)
        
        return f
    
    dzeta = a8*f0(a9*(a10-gsi))+a12*f0(a13*(a14-gsi))+a16*f0(a17*(a18-gsi))+(a13*math.exp(a13*(a14-gsi))*(a11+a12*gsi))*math.pow(f0(a13*(a14-gsi)),2)+(a17*math.exp(a17*(a18-gsi))*(a15+a16*gsi))*math.pow(f0(a17*(a18-gsi)),2)+(a9*math.exp(a9*(a10-gsi))*(a7+a8*gsi))*math.pow(f0(a9*(a10-gsi)),2)+(a2+3.0*a3*math.pow(gsi,2))*f0(a5*(gsi-a6))/(1+a4*gsi)-(a4*(a1+a2*gsi+a3*math.pow(gsi,3)))*f0(a5*(gsi-a6))/math.pow((1.0+a4*gsi),2)-(a5*math.exp(a5*(gsi-a6))*(a1+a2*gsi+a3*math.pow(gsi,3)))*math.pow(f0(a5*(gsi-a6)),2)/(1.0+a4*gsi)
    
    return dzeta

def fps_d2zeta(gsi):
    
    a1 = 6.22
    a2 = 6.121
    a3 = 0.006004
    a4 = 0.16345
    a5 = 6.50
    a6 = 11.8440
    a7 = 17.24
    a8 = 1.065
    a9 = 6.54
    a10 = 11.8421
    a11 = -22.003
    a12 = 1.5552
    a13 = 9.3
    a14 = 14.19
    a15 = 23.73
    a16 = -1.508
    a17 = 1.79
    a18 = 15.13
    
    # defining the function f0
    # returns the value of f0 at x
    def f0(x):
        
        f = 1.0/(math.exp(x)+1)
        
        return f
    
    ddzeta = (2.0*a8*a9*math.exp(a9*(a10-gsi)))*math.pow(f0(a9*(a10-gsi)),2)+(2.0*a12*a13*math.exp(a13*(a14-gsi)))*math.pow(f0(a13*(a14-gsi)),2)+(2.0*a16*a17*math.exp(a17*(a18-gsi)))*math.pow(f0(a17*(a18-gsi)),2)+(2.0*math.pow(a13,2)*math.exp(2*a13*(a14-gsi))*(a11+a12*gsi))*math.pow(f0(a13*(a14-gsi)),3)-(math.pow(a13,2)*math.exp(a13*(a14-gsi))*(a11+a12*gsi))*math.pow(f0(a13*(a14-gsi)),2)+(2.0*math.pow(a17,2)*math.exp(2.0*a17*(a18-gsi))*(a15+a16*gsi))*math.pow(f0(a17*(a18-gsi)),3)-(math.pow(a17,2)*math.exp(a17*(a18-gsi))*(a15+a16*gsi))*math.pow(f0(a17*(a18-gsi)),2)+(6.0*a3*gsi)*f0(a5*(-a6+gsi))/(1.0+a4*gsi)+(2.0*math.pow(a9,2)*math.exp(2.0*a9*(a10-gsi))*(a7+a8*gsi))*math.pow(f0(a9*(a10-gsi)),3)-(math.pow(a9,2)*math.exp(a9*(a10-gsi))*(a7+a8*gsi))*math.pow(f0(a9*(a10-gsi)),2)-(2.0*a4*(a2+3.0*a3*math.pow(gsi,2)))*f0(a5*(-a6+gsi))/math.pow((1.0+a4*gsi),2)-(2.0*a5*math.exp(a5*(-a6+gsi))*(a2+3.0*a3*math.pow(gsi,2)))*math.pow(f0(a5*(-a6+gsi)),2)/(1+a4*gsi)+(2.0*math.pow(a4,2)*(a1+a2*gsi+a3*math.pow(gsi,3)))*f0(a5*(-a6+gsi))/math.pow((1.0+a4*gsi),3)+(2.0*a4*a5*math.exp(a5*(-a6+gsi))*(a1+a2*gsi+a3*math.pow(gsi,3)))*math.pow(f0(a5*(-a6+gsi)),2)/math.pow((1.0+a4*gsi),2)+(2.0*math.pow(a5,2)*math.exp(2.0*a5*(-a6+gsi))*(a1+a2*gsi+a3*math.pow(gsi,3)))*math.pow(f0(a5*(-a6+gsi)),3)/(1.0+a4*gsi)-(math.pow(a5,2)*math.exp(a5*(-a6+gsi))*(a1+a2*gsi+a3*math.pow(gsi,3)))*math.pow(f0(a5*(-a6+gsi)),2)/(1.0+a4*gsi)
    
    return ddzeta
    

print("These are a few, of my most favorite things")
