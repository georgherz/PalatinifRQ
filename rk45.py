import numpy as np
import pdb



def rk4_step(fx, yn, t, h, tol, rho):#, EOS, conv): 
    #makes one rkf45 step
    # fx is the function which should be integrated (has to be in the form dx = ...), also a vector of functions can be passed on.
    # yn is the vector of initial conditions. length of yn vector and fx vector need to be the same
    # t is the parameter according to which I integrate my function (i.e. I have dx/dt --> I integrate according to t)
    # h is my stepsize
    # tol is my error tolerance
    # rho is my vector of parameters I have to pass on
    
    n = len(yn) # saves the length of my initial condition vector
    
    k1 = np.zeros(n) # creats a vector of length n filled with zeros
    k2 = np.zeros(n) # see above
    k3 = np.zeros(n) # see above
    k4 = np.zeros(n) # see above
    k5 = np.zeros(n) # see above
    k6 = np.zeros(n) # see above
    yn1 = np.zeros(n) # see above (will be the value at n+1 for the rk4 step)
    zn1 = np.zeros(n) # see above (will be the value at n+1 for the rk5 step)
    
    #pdb.set_trace()
    
    # Calculate the 6 trial steps necessary for the rkf45 algorythm
    k1 = h*fx(t,yn,rho) # the functions stored in the vector fx are calculated at the point t and yn and then saved at k1
    k2 = h*fx(t+1.0/4.0*h,yn+1.0/4.0*k1,rho) # same as above with t+1/4*h and yn+1/4*k1
    k3 = h*fx(t+3.0/8.0*h,yn+3.0/32.0*k1+9.0/32.0*k2,rho) # same as above with t+3/8*h and yn+3/32*k1+9/32*k2
    k4 = h*fx(t+12.0/13.0*h,yn+1932.0/2197.0*k1-7200.0/2197.0*k2+7296.0/2197.0*k3,rho) # calculating k4 for runge-kutta-fehlberg method
    k5 = h*fx(t+h, yn+439.0/216.0*k1-8*k2+3680.0/513.0*k3-845.0/4104.0*k4,rho) # k5 of rkf45 method
    k6 = h*fx(t+1.0/2.0*h, yn-8.0/27.0*k1+2*k2-3544.0/2565.0*k3+1859.0/4104.0*k4-11.0/40.0*k5,rho) # k6 of rkf45 method
    
    #pdb.set_trace()
    
    # use an RK4 step to approximate the the solution    
    yn1 = yn+25.0/216.0*k1+1408.0/2565.0*k3+2197.0/4101.0*k4-1.0/5.0*k5
    # use an RK5 step to determine a better value for the solution
    zn1 = yn+16.0/135.0*k1+6656.0/12825.0*k3+28561.0/56430.0*k4-9.0/50.0*k5+2.0/55.0*k6
    
    # determine the optimal stepsize
    H_ideal = h*np.power(np.absolute(yn1-zn1)/tol,-0.20)
    h_ideal = H_ideal.min() # graps the minimal stepsize of all the possible stepsizes

    
    if h_ideal < h:
        return np.array([yn,h_ideal*0.75,t])
    else:
        return np.array([yn1, h_ideal, t+h])
                         
    
    

print("Stimmungshochhalter")
