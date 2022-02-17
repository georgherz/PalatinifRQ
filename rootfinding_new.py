import numpy as np
import math
import scipy
from scipy import optimize
#import pdb


def rootfinding(fx, y, init_guess):
    # Calculates the inverse of the function fx by finding the roots of h(x) = fx - y
    # The roots of h(x) are the values x for which fx gives y, i.e. fx(root) = y
    # Input: fx ...  the function I want to numerically invert, y ... the value of fx at x, init_guess ... the initial guess for the root
    # returns the root of h(x) which is the inverse of fx
    # the newton method from scipy.optimize is used to calculate the roots of h(x)

    def h(x):   # defining the function h(x)
        h = fx(x)-y
        return h
    
    out = optimize.newton(h,init_guess) # calculates the roots of h via the newton method  
    return out

print("These are a few, of my most favorite things")
