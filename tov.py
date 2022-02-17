import numpy as np
import math

def TOV(t, yn, rho):
    
    # the rhs of all functions I need to integrate in order to solve the TOV equation are given here
    
    r = t # defines my integration parameter as r (so that I am less confused when writing the equations
    dm = 4*math.pi*rho*pow(r,2) # equation for mass according to Wald 6.2.10
    dp = -(yn[0]+rho)*(yn[1]+4*math.pi*pow(r,3)*yn[0])/(r*(r-2*yn[1])) # TOV equation according to Wald 6.2.19
    return np.array([dp, dm]) # returns my functions in an array
