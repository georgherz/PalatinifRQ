import numpy as np
import math
import pdb
import matplotlib.pylab as plt


def rk4_step(fx, yn, t, h, rho): 
    #makes one rk4 step
    # fx is the function which should be integrated (has to be in the form dx = ...), also a vector of functions can be passed on.
    # yn is the vector of initial conditions. length of yn vector and fx vector need to be the same
    # t is the parameter according to which I integrate my function (i.e. I have dx/dt --> I integrate according to t)
    # h is my stepsize
    
    n = len(yn) # saves the length of my initial condition vector
    
    k1 = np.zeros(n) # creats a vector of length n filled with zeros
    k2 = np.zeros(n) # see above
    k3 = np.zeros(n) # see above
    k4 = np.zeros(n) # see above
    yn1 = np.zeros(n) # see above
    #pdb.set_trace()
    k1 = fx(t,yn, rho) # the functions stored in the vector fx are calculated at the point t and yn and then saved at k1
    k2 = fx(t+0.5*h,yn+0.5*h*k1, rho) # same as above with t+0.5*h and yn+0.5*h*k1
    k3 = fx(t+0.5*h,yn+0.5*h*k2, rho) # same as above with t+0.5*h and yn+0.5*h*k2
    k4 = fx(t+h,yn+h*k3, rho) # same as above with t+h and yn+h*k3
    
    yn1 = yn + h/6.0*(k1+2*(k2+k3)+k4) # makes the rk4 step
    
    #print(k1)
    
    return yn1 # new values for all the functions I am integrating are returned
