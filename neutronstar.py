import numpy as np
import math
import rootfinding_new as rf
import sly_ply_fps_new as spf
import matplotlib.pylab as plt
from scipy.interpolate import interp1d
import rk45 as rkf
import rk4
import tov
import time
import modified_tov_fR as mod_tov_fR
import modified_tov_fRQ as mod_tov_fRQ
import modified_tov_tab_fRQ as mod_tov_tab_fRQ
import modified_tov_tab_fR as mod_tov_tab_fR
import pdb

### Conversion factors cgs --> geometrized
rho_conv = 7.426*math.pow(10, -29)
p_conv = 8.262*math.pow(10, -50)
MeV_conv = 1.602*math.pow(10,33) # converts our rho and p from MeV/fm^3 to baryes

Eos = [spf.fps_zeta, "HS_DD2.dat"] # my equations of state
h = 100.0 # stepsize
alfa = 1.0*np.power(10.0, 9) # the alpha I am using for f(R)
beta = 1.0*np.power(10.0, 9)
gsi_c = 14.5 # set the value for the lowest gsi we are considering
uconv = [p_conv, rho_conv] # put the conversion factors into an array

#### GR #####
RK4_results = []
myfile = open("GR.dat", "w")
radius_HSDD2_GR = []
gsi_HSDD2_GR = []
mass_HSDD2_GR = []
radius_FPS_GR = []
gsi_FPS_GR = []
mass_FPS_GR = []
for fkt in Eos:
    gsi_c=14.5
    if fkt == "HS_DD2.dat":
        ### Make my interpolation functions
        #### Read data
        eos = open(fkt, "r") # Opens the data file
        eos_data = eos.readlines() # reads in the lines of the data file and saves them in eos_data --> all the data is in str format
            
        #### Converts the data from str into float (if necessary)
        eos_useful = [] # empty array to store the data in
        for idx in range(len(eos_data)): # loops over the full array
            eos_useful.append(eos_data[idx].split()) # splits the string in eos_data whenever there is a space and then appends it to list --> we get a list where each entry has 5 entries ([gsi], [m], [r], [EOS], [alpha])
            eos_useful[idx][0]=float(eos_useful[idx][0]) # converts from str to float
            eos_useful[idx][1]=float(eos_useful[idx][1]) # same as above
            eos_useful[idx][2]=float(eos_useful[idx][2]) # same as above
            eos_useful[idx][3]=float(eos_useful[idx][3]) # same as above
            
        #### Arrange all values in a way that temperature, pressure, etc all are in one array
        EosData = np.column_stack(eos_useful[0:]) # since we calculated the derivatives without any unit conversions I have to convert everything into geometrized units (see next lines)
        rho = EosData[1]*MeV_conv*p_conv # calculate the energydensity in geometrized units. Energydensity has the same conversion factor as pressure
        pressure = EosData[0]*MeV_conv*p_conv # calculate the pressure in geometrized units
        drho = EosData[2] # conversion factors for depsilon cancel our
        ddrho = EosData[3]/(MeV_conv*p_conv) # calculate ddepsilon in geometrized units
            
        pinterpol = interp1d(rho, pressure) # interpolates my pressure for a given density, i.e. putting q1 into pinterpol gives me the pressure associated with the density q1 by interpolation
        qinterpol = interp1d(pressure, rho) # Interpolates my density for a given pressure, i.e. qinterpol(p1) gives me the density associated with p1 by interpolation
        dqinterpol = interp1d(pressure, drho) # Interpolation for drho
        ddqinterpol = interp1d(pressure, ddrho) # Interpolation for ddrho
    
    if fkt == "HS_DD2.dat":
        gsi_final = 15.7
    else:
        gsi_final = 16
    
    while gsi_c < gsi_final:
        m0 = 0 # initialize the mass (mass at the core of the neutron star is zero)
        if fkt == "HS_DD2.dat":
            p_c = pinterpol(math.pow(10, gsi_c)*uconv[1]) # calculate the central pressure via interpolation
            print(str(gsi_c)+" "+str(p_c)+" "+str(fkt))
        else:
            p_c = math.pow(10,fkt(gsi_c))*uconv[0] # calculate the central pressure
            print(str(gsi_c)+" "+str(p_c)+" "+str(fkt))
        t = 0 # set the integration parameter (will be the radius later) to zero
        if fkt == "HS_DD2.dat":
            p_stop = pressure[0]
        else:
            p_stop = p_c*math.pow(10,-12) # calculate the pressure when we stop the integration, we stop when the pressure is 12 orders of magnitude lower than the initial pressure
        y0 = np.array([p_c, m0]) # array with initial values
        r = (t+0.0001)*h # start the integration process a little bit off the center in order to avoid the numerical singularity at r=0
        #parameters for fR
        starttime = time.time() # calculate the starting time and save it
        # as long as the pressure is bigger than our final pressure, perform the integration
        while y0[0] > p_stop:
            if fkt == "HS_DD2.dat":
                try: 
                    rho = qinterpol(y0[0])
                    out = rk4.rk4_step(tov.TOV, y0, r, h, rho) # makes the RKF45 step and saves the new pressure, mass, ideal stepsize and radious in the variable out
                    y0[0] = out[0]
                    y0[1] = out[1]
                    r = r+h
                except ValueError:
                    break
            else:
                gsi = rf.rootfinding(fkt,np.log10(y0[0]/uconv[0]),0) # find the gsi associated with the pressure numerically
                rho = math.pow(10,gsi)*uconv[1] #calculate rho (I go into the TOV eq with it)
                out = rk4.rk4_step(tov.TOV, y0, r, h, rho) # makes the RK4 step and saves the new pressure, mass, ideal stepsize and radious in the variable out
                y0[0] = out[0]
                y0[1] = out[1]
                r = r+h

        if fkt == "HS_DD2.dat":
            radius_HSDD2_GR.append(r)
            gsi_HSDD2_GR.append(gsi_c)
            mass_HSDD2_GR.append(y0[1])
        else:
            radius_FPS_GR.append(r)
            gsi_FPS_GR.append(gsi_c)
            mass_FPS_GR.append(y0[1])
        endtime = time.time() # calculates the time at the end of integration
        # print elapsed time to monitor integration speed
        print("Elapsed time in seconds for stepsize h="+str(h))
        print(endtime-starttime)
        fulltime = endtime-starttime
        # print the final values
        print("Final Values")
        print(r)
        print(y0[1])
            
        # depending on which EOS I am using, append the results to different arrays
        if fkt == "HS_DD2.dat":
            RK4_results.append([y0[1]/(uconv[1]*1.99*math.pow(10,33)), r/100000.0, fulltime, h, "HSDD2", alfa, "fR"])
            myfile.write(str(gsi_c)+" "+str(y0[1]/(uconv[1]*1.99*math.pow(10,33)))+" "+str(r/100000.0)+" "+str(fulltime)+" "+str(h)+" HSDD2 GR "+"\n")
        if fkt == spf.fps_zeta:
            RK4_results.append([y0[1]/(uconv[1]*1.99*math.pow(10,33)), r/100000.0, fulltime, h, "FPS", alfa, "fR"])
            myfile.write(str(gsi_c)+" "+str(y0[1]/(uconv[1]*1.99*math.pow(10,33)))+" "+str(r/100000.0)+" "+str(fulltime)+" "+str(h)+" FPS GR "+"\n")
        
        gsi_c = gsi_c+0.1

myfile.close()

#### fR #####
RK4_results = []
myfile = open("fR.dat", "w")
radius_HSDD2_fR = []
gsi_HSDD2_fR = []
mass_HSDD2_fR = []
radius_FPS_fR = []
gsi_FPS_fR = []
mass_FPS_fR = []
for fkt in Eos:
    gsi_c = 14.5
    if fkt == "HS_DD2.dat":
            
        ### Make my interpolation functions
        #### Read data
        eos = open(fkt, "r") # Opens the data file
        eos_data = eos.readlines() # reads in the lines of the data file and saves them in eos_data --> all the data is in str format
            
        #### Converts the data from str into float (if necessary)
        eos_useful = [] # empty array to store the data in
        for idx in range(len(eos_data)): # loops over the full array
            eos_useful.append(eos_data[idx].split()) # splits the string in eos_data whenever there is a space and then appends it to list --> we get a list where each entry has 5 entries ([gsi], [m], [r], [EOS], [alpha])
            eos_useful[idx][0]=float(eos_useful[idx][0]) # converts from str to float
            eos_useful[idx][1]=float(eos_useful[idx][1]) # same as above
            eos_useful[idx][2]=float(eos_useful[idx][2]) # same as above
            eos_useful[idx][3]=float(eos_useful[idx][3]) # same as above
            
        #### Arrange all values in a way that temperature, pressure, etc all are in one array
        EosData = np.column_stack(eos_useful[0:]) # since we calculated the derivatives without any unit conversions I have to convert everything into geometrized units (see next lines)
        rho = EosData[1]*MeV_conv*p_conv # calculate the energydensity in geometrized units. Energydensity has the same conversion factor as pressure
        pressure = EosData[0]*MeV_conv*p_conv # calculate the pressure in geometrized units
        drho = EosData[2] # conversion factors for depsilon cancel our
        ddrho = EosData[3]/(MeV_conv*p_conv) # calculate ddepsilon in geometrized units
            
        pinterpol = interp1d(rho, pressure) # interpolates my pressure for a given density, i.e. putting q1 into pinterpol gives me the pressure associated with the density q1 by interpolation
        qinterpol = interp1d(pressure, rho) # Interpolates my density for a given pressure, i.e. qinterpol(p1) gives me the density associated with p1 by interpolation
        dqinterpol = interp1d(pressure, drho) # Interpolation for drho
        ddqinterpol = interp1d(pressure, ddrho) # Interpolation for ddrho
    if fkt == "HS_DD2.dat":
        gsi_final = 15.7
    else:
        gsi_final = 16
        
    while gsi_c < gsi_final:
        m0 = 0 # initialize the mass (mass at the core of the neutron star is zero)
        if fkt == "HS_DD2.dat":
            p_c = pinterpol(math.pow(10, gsi_c)*uconv[1]) # calculate the central pressure via interpolation
            print(str(gsi_c)+" "+str(p_c)+" "+str(fkt))
        else:
            p_c = math.pow(10,fkt(gsi_c))*uconv[0] # calculate the central pressure
            print(str(gsi_c)+" "+str(p_c)+" "+str(fkt))
        t = 0 # set the integration parameter (will be the radius later) to zero
        if fkt == "HS_DD2.dat":
            p_stop = pressure[0]
        else:
            p_stop = p_c*math.pow(10,-12) # calculate the pressure when we stop the integration, we stop when the pressure is 12 orders of magnitude lower than the initial pressure
        y0 = np.array([p_c, m0]) # array with initial values
        r = (t+0.0001)*h # start the integration process a little bit off the center in order to avoid the numerical singularity at r=0
        #parameters for fR
        if fkt == "HS_DD2.dat":
            para = [0.0, 8*math.pi, alfa, dqinterpol, ddqinterpol] # rho (is only set to zero for initializing, in the next loop it will be set to its value and adapted everytime), kappa, alpha1, dEOS, d2EOS
        if fkt == spf.fps_zeta:
            para = [0.0, 8*math.pi, alfa, spf.fps_dzeta, spf.fps_d2zeta] # rho (is only set to zero for initializing, in the next loop it will be set to its value and adapted everytime), kappa, alpha1, dEOS, d2EOS
        starttime = time.time() # calculate the starting time and save it
        # as long as the pressure is bigger than our final pressure, perform the integration
        while y0[0] > p_stop:
            if fkt == "HS_DD2.dat":
                try: 
                    rho = qinterpol(y0[0])
                    para[0]= rho
                    out = rk4.rk4_step(mod_tov_tab_fR.modified_TOV, y0, r, h, para) # makes the RKF45 step and saves the new pressure, mass, ideal stepsize and radious in the variable out
                    p_old = y0[0]
                    m_old = y0[1]
                    y0[0] = out[0]
                    y0[1] = out[1]
                    r = r+h
                except ValueError:
                    break
            else:
                gsi = rf.rootfinding(fkt,np.log10(y0[0]/uconv[0]),0) # find the gsi associated with the pressure numerically
                rho = math.pow(10,gsi)*uconv[1] #calculate rho (I go into the TOV eq with it)
                para[0]= rho
                out = rk4.rk4_step(mod_tov_fR.modified_TOV, y0, r, h, para) # makes the RK4 step and saves the new pressure, mass, ideal stepsize and radious in the variable out
                y0[0] = out[0]
                y0[1] = out[1]
                r = r+h
        if fkt == "HS_DD2.dat":
            radius_HSDD2_fR.append(r)
            gsi_HSDD2_fR.append(gsi_c)
            mass_HSDD2_fR.append(y0[1])
        else:
            radius_FPS_fR.append(r)
            gsi_FPS_fR.append(gsi_c)
            mass_FPS_fR.append(y0[1])
        endtime = time.time() # calculates the time at the end of integration
        # print elapsed time to monitor integration speed
        print("Elapsed time in seconds for stepsize h="+str(h))
        print(endtime-starttime)
        fulltime = endtime-starttime
        # print the final values
        print("Final Values")
        print(r)
        print(y0[1])
            
        # depending on which EOS I am using, append the results to different arrays
        if fkt == "HS_DD2.dat":
            RK4_results.append([y0[1]/(uconv[1]*1.99*math.pow(10,33)), r/100000.0, fulltime, h, "HSDD2", alfa, "fR"])
            myfile.write(str(gsi_c)+" "+str(y0[1]/(uconv[1]*1.99*math.pow(10,33)))+" "+str(r/100000.0)+" "+str(fulltime)+" "+str(h)+" HSDD2 fR "+str(alfa)+"\n")
        if fkt == spf.fps_zeta:
            RK4_results.append([y0[1]/(uconv[1]*1.99*math.pow(10,33)), r/100000.0, fulltime, h, "FPS", alfa, "fR"])
            myfile.write(str(gsi_c)+" "+str(y0[1]/(uconv[1]*1.99*math.pow(10,33)))+" "+str(r/100000.0)+" "+str(fulltime)+" "+str(h)+" FPS fR "+str(alfa)+"\n")
            
        gsi_c = gsi_c+0.1

myfile.close()

### fRQ ###
RK4_results = []
myfile = open("fRQ.dat", "w")
radius_HSDD2_fRQ = []
gsi_HSDD2_fRQ = []
mass_HSDD2_fRQ = []
radius_FPS_fRQ = []
gsi_FPS_fRQ = []
mass_FPS_fRQ = []
for fkt in Eos:
    gsi_c = 14.5
    if fkt == "HS_DD2.dat":
        
        ### Make my interpolation functions
        #### Read data
        eos = open(fkt, "r") # Opens the data file
        eos_data = eos.readlines() # reads in the lines of the data file and saves them in eos_data --> all the data is in str format
        
        #### Converts the data from str into float (if necessary)
        eos_useful = [] # empty array to store the data in
        for idx in range(len(eos_data)): # loops over the full array
            eos_useful.append(eos_data[idx].split()) # splits the string in eos_data whenever there is a space and then appends it to list --> we get a list where each entry has 5 entries ([gsi], [m], [r], [EOS], [alpha])
            eos_useful[idx][0]=float(eos_useful[idx][0]) # converts from str to float
            eos_useful[idx][1]=float(eos_useful[idx][1]) # same as above
            eos_useful[idx][2]=float(eos_useful[idx][2]) # same as above
            eos_useful[idx][3]=float(eos_useful[idx][3]) # same as above
            
        #### Arrange all values in a way that temperature, pressure, etc all are in one array
        EosData = np.column_stack(eos_useful[0:]) # since we calculated the derivatives without any unit conversions I have to convert everything into geometrized units (see next lines)
        rho = EosData[1]*MeV_conv*p_conv # calculate the energydensity in geometrized units. Energydensity has the same conversion factor as pressure
        pressure = EosData[0]*MeV_conv*p_conv # calculate the pressure in geometrized units
        drho = EosData[2] # conversion factors for depsilon cancel our
        ddrho = EosData[3]/(MeV_conv*p_conv) # calculate ddepsilon in geometrized units
        
        pinterpol = interp1d(rho, pressure) # interpolates my pressure for a given density, i.e. putting q1 into pinterpol gives me the pressure associated with the density q1 by interpolation
        qinterpol = interp1d(pressure, rho) # Interpolates my density for a given pressure, i.e. qinterpol(p1) gives me the density associated with p1 by interpolation
        dqinterpol = interp1d(pressure, drho) # Interpolation for drho
        ddqinterpol = interp1d(pressure, ddrho) # Interpolation for ddrho
    
    if fkt == "HS_DD2.dat":
        gsi_final = 15.7
    else:
        gsi_final = 16
            
    while gsi_c < gsi_final:
        m0 = 0 # initialize the mass (mass at the core of the neutron star is zero)
        if fkt == "HS_DD2.dat":
            p_c = pinterpol(math.pow(10, gsi_c)*uconv[1]) # calculate the central pressure via interpolation
            print(str(gsi_c)+" "+str(p_c)+" "+str(fkt))
        else:
            p_c = math.pow(10,fkt(gsi_c))*uconv[0] # calculate the central pressure
            print(str(gsi_c)+" "+str(p_c)+" "+str(fkt))
        t = 0 # set the integration parameter (will be the radius later) to zero
        if fkt == "HS_DD2.dat":
            p_stop = pressure[0]
        else:
            p_stop = p_c*math.pow(10,-12) # calculate the pressure when we stop the integration, we stop when the pressure is 12 orders of magnitude lower than the initial pressure
        y0 = np.array([p_c, m0]) # array with initial values
        r = (t+0.0001)*h # start the integration process a little bit off the center in order to avoid the numerical singularity at r=0
        #parameters for fR
        if fkt == "HS_DD2.dat":
            para = [0.0, 8*math.pi, alfa, beta, dqinterpol, ddqinterpol] # rho (is only set to zero for initializing, in the next loop it will be set to its value and adapted everytime), kappa, alpha1, dEOS, d2EOS
        if fkt == spf.fps_zeta:
            para = [0.0, 8*math.pi, alfa, beta, spf.fps_dzeta, spf.fps_d2zeta] # rho (is only set to zero for initializing, in the next loop it will be set to its value and adapted everytime), kappa, alpha1, beta1, dEOS, d2EOS
        starttime = time.time() # calculate the starting time and save it
        # as long as the pressure is bigger than our final pressure, perform the integration
        while y0[0] > p_stop:
            if fkt == "HS_DD2.dat":
                try: 
                    rho = qinterpol(y0[0])
                    para[0]= rho
                    out = rk4.rk4_step(mod_tov_tab_fRQ.modified_TOV, y0, r, h, para) # makes the RKF45 step and saves the new pressure, mass, ideal stepsize and radious in the variable out
                    p_old = y0[0]
                    m_old = y0[1]
                    y0[0] = out[0]
                    y0[1] = out[1]
                    r = r+h
                except ValueError:
                    break
            else:
                gsi = rf.rootfinding(fkt,np.log10(y0[0]/uconv[0]),0) # find the gsi associated with the pressure numerically
                rho = math.pow(10,gsi)*uconv[1] #calculate rho (I go into the TOV eq with it)
                para[0]= rho
                out = rk4.rk4_step(mod_tov_fRQ.modified_TOV, y0, r, h, para) # makes the RK4 step and saves the new pressure, mass, ideal stepsize and radious in the variable out
                y0[0] = out[0]
                y0[1] = out[1]
                r = r+h
        if fkt == "HS_DD2.dat":
            radius_HSDD2_fRQ.append(r)
            gsi_HSDD2_fRQ.append(gsi_c)
            mass_HSDD2_fRQ.append(y0[1])
        else:
            radius_FPS_fRQ.append(r)
            gsi_FPS_fRQ.append(gsi_c)
            mass_FPS_fRQ.append(y0[1])
        endtime = time.time() # calculates the time at the end of integration
        # print elapsed time to monitor integration speed
        print("Elapsed time in seconds for stepsize h="+str(h))
        print(endtime-starttime)
        fulltime = endtime-starttime
        # print the final values
        print("Final Values")
        print(r)
        print(y0[1])
        
        # depending on which EOS I am using, append the results to different arrays
        if fkt == "HS_DD2.dat":
            RK4_results.append([y0[1]/(uconv[1]*1.99*math.pow(10,33)), r/100000.0, fulltime, h, "HSDD2", alfa, beta, "fRQ"])
            myfile.write(str(gsi_c)+" "+str(y0[1]/(uconv[1]*1.99*math.pow(10,33)))+" "+str(r/100000.0)+" "+str(fulltime)+" "+str(h)+" HSDD2 fRQ "+str(alfa)+" "+str(beta)+"\n")
        if fkt == spf.fps_zeta:
            RK4_results.append([y0[1]/(uconv[1]*1.99*math.pow(10,33)), r/100000.0, fulltime, h, "FPS", alfa, beta, "fRQ"])
            myfile.write(str(gsi_c)+" "+str(y0[1]/(uconv[1]*1.99*math.pow(10,33)))+" "+str(r/100000.0)+" "+str(fulltime)+" "+str(h)+" FPS fRQ "+str(alfa)+" "+str(beta)+"\n")
            
        gsi_c = gsi_c+0.1

myfile.close()

# mass-radius plot
plt.figure(dpi=300)
plt.plot(np.array(radius_HSDD2_GR)/100000, np.array(mass_HSDD2_GR)/(uconv[1]*1.99*math.pow(10,33)), label="HS(DD2) GR", linewidth = 0.5, color="black")
plt.plot(np.array(radius_HSDD2_fR)/100000, np.array(mass_HSDD2_fR)/(uconv[1]*1.99*math.pow(10,33)), label="HS(DD2) f(R)", linewidth = 0.5, color="blue")
plt.plot(np.array(radius_HSDD2_fRQ)/100000, np.array(mass_HSDD2_fRQ)/(uconv[1]*1.99*math.pow(10,33)), label="HS(DD2) f(R,Q)", linewidth = 0.5, color="orange")
plt.plot(np.array(radius_FPS_GR)/100000, np.array(mass_FPS_GR)/(uconv[1]*1.99*math.pow(10,33)), label="FPS GR", linewidth = 0.5, color="black", linestyle="--")
plt.plot(np.array(radius_FPS_fR)/100000, np.array(mass_FPS_fR)/(uconv[1]*1.99*math.pow(10,33)), label="FPS f(R)", linewidth = 0.5, color="blue", linestyle="--")
plt.plot(np.array(radius_FPS_fRQ)/100000, np.array(mass_FPS_fRQ)/(uconv[1]*1.99*math.pow(10,33)), label="FPS f(R,Q)", linewidth = 0.5, color="orange", linestyle="--")
plt.legend(prop={'size': 8}, loc="upper right")
plt.xlabel(r'$R\ [km]$')
plt.ylabel(r'M $[m_{\odot}]$')
plt.xlim(8,20)
plt.ylim(0,4)
plt.savefig("m_vs_r_exampleCode.png")
plt.clf()


# gsi_c-mass plot
plt.figure(dpi=300)
plt.plot(np.array(gsi_HSDD2_GR), np.array(mass_HSDD2_GR)/(uconv[1]*1.99*math.pow(10,33)), label="GR HS(DD2)", linewidth = 0.5, color = "black")
plt.plot(np.array(gsi_HSDD2_fR), np.array(mass_HSDD2_fR)/(uconv[1]*1.99*math.pow(10,33)), label="HS(DD2) f(R)", linewidth = 0.5, color="blue")
plt.plot(np.array(gsi_HSDD2_fRQ), np.array(mass_HSDD2_fRQ)/(uconv[1]*1.99*math.pow(10,33)), label="HS(DD2) f(R,Q)", linewidth = 0.5, color="orange")
plt.plot(np.array(gsi_FPS_GR), np.array(mass_FPS_GR)/(uconv[1]*1.99*math.pow(10,33)), label="FPS GR", linewidth = 0.5, color="black", linestyle="--")
plt.plot(np.array(gsi_FPS_fR), np.array(mass_FPS_fR)/(uconv[1]*1.99*math.pow(10,33)), label="FPS f(R)", linewidth = 0.5, color="blue", linestyle="--")
plt.plot(np.array(gsi_FPS_fRQ), np.array(mass_FPS_fRQ)/(uconv[1]*1.99*math.pow(10,33)), label="FPS f(R,Q)", linewidth = 0.5, color="orange", linestyle="--")
plt.legend(prop={'size': 8}, loc="upper right")
plt.xlabel(r'$log_{10} \rho_c\ [\frac{g}{cm^3}]$')
plt.ylabel(r'M $[m_{\odot}]$')
plt.xlim(14,16.5)
plt.ylim(0,4)
plt.savefig("m_vs_rho_c_exampleCode.png")
plt.clf()
