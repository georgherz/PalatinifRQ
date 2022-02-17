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
import modified_tov_tab_fRQ as mod_tov_fRQ
import pdb
    

### Set initial conditions ###

#myfunctarray = [spf.ply_zeta, spf.sly_zeta, spf.fps_zeta] # my equations of state
myfunctarray = ["HS_DD2.dat"] # my equations of state
rho_conv = 7.426*math.pow(10, -29)
p_conv = 8.262*math.pow(10, -50)
MeV_conv = 1.602*math.pow(10,33) # converts our rho and p from MeV/fm^3 to baryes

tol = math.pow(10, 1) # sets the error tolerance between two rkf45 steps


# make a loop over all alpha1 (couping constants for the R^2 term)
alphaarray = [-2.0*np.power(10.0, 9), -4.0*np.power(10.0, 9), -1.0*np.power(10.0, 8), 4.0*np.power(10.0, 9), 3.0*np.power(10.0, 9), 1.0*np.power(10.0, 9), 2.0*np.power(10.0, 9), -1.0*np.power(10.0, 9), 1.0*np.power(10.0, 8), 2.0*np.power(10.0, -8)]
betaarray = [1.0*np.power(10.0,9), 2.0*np.power(10.0,9), 3.0*np.power(10.0,9), 4.0*np.power(10.0,9), 5.0*np.power(10.0,9), 1.0*np.power(10.0, -8)]
alphaarray = [1.0*np.power(10.0,-8), 1.0*np.power(10.0,9)]
betaarray = [1.0*np.power(10.0,9), 4.0*np.power(10,9)]

h = 100.0
fulltime1 = time.time()

# initialize the arrays where the full data is saved
full_HSDD2 = []
full_sly = []
full_fps = []
fullfR_HSDD2 = []
fullfR_sly = []
fullfR_fps = []

myfile = open("stardata.dat", "w")

for alfa in alphaarray: # loop over all alpha
    
    for beta in betaarray: # loop over all beta
    
        # initialize the arrays where the final values after each run are saved
        printimass_HSDD2 = [] # saves the mass
        printimass_sly = [] # same as above
        printimass_fps = [] # same as above
        printifR_HSDD2 = []
        printifR_sly = []
        printifR_fps = []
        printiradius_HSDD2 = [] # saves the radius
        printiradius_sly = [] # saves the radius
        printiradius_fps = [] # saves the radius
        printifRradius_HSDD2 = [] # saves the radii for fR vs r, i.e. the full radius from center to surface
        printifRradius_sly = [] # same as above
        printifRradius_fps = [] # same as above
        printi_rho_HSDD2 = []
        printi_rho_sly = []
        printi_rho_fps = []


        for fkt in myfunctarray: # make a loop over all possible equations of state
            
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
            
            gsi_c = 14.4 # set the value for the lowest gsi we are considering
            gsi_final = 15.6
            printigsi = [] # initialize the array where we save all the gsi values we used in
            uconv = [p_conv, rho_conv] # put the conversion factors into an array
    
            while gsi_c < gsi_final: # loop over all values of gsi until gsi = 16
                
                m0 = 0 # initialize the mass (mass at the core of the neutron star is zero)
                #p_c = math.pow(10,fkt(gsi_c))*uconv[0] # calculate the central pressure
                p_c = pinterpol(math.pow(10, gsi_c)*uconv[1])
                t = 0 # set the integration parameter (will be the radius later) to zero
                #p_stop = p_c*math.pow(10,-12) # calculate the pressure when we stop the integration, we stop when the pressure is 12 orders of magnitude lower than the initial pressure
                p_stop = pressure[0]
                y0 = np.array([p_c, m0]) # array with initial values
    
                # print my starting values so that I can monitor progress
                print("Start Values")
                print("Central density: "+str(gsi_c))
                print("Alpha: "+str('{:.1e}'.format(alfa)))
                print("Beta: "+str('{:.1e}'.format(beta)))
                print("Eos: "+str(fkt))
                
                r = (t+0.0001)*h # start the integration process a little bit off the center in order to avoid the numerical singularity at r=0
            
                starttime = time.time() # calculate the starting time and save it
            
                #para = [0.0, 8.0*math.pi*6.67*np.power(10.0, -11), 2.0*np.power(10.0, -8), 3.0*np.power(10.0, 9), 2.0*np.power(10.0, -8), spf.ply_dzeta, spf.ply_d2zeta] # rho (is only set to zero for initializing, in the next loop it will be set to its value and adapted everytime), kappa, alpha1, beta1, fQ, dEOS, d2EOS
            
                #parameters for fRQ
                if fkt == "HS_DD2.dat":
                    para = [0.0, 8*math.pi, alfa, beta, dqinterpol, ddqinterpol] # rho (is only set to zero for initializing, in the next loop it will be set to its value and adapted everytime), kappa, alpha1, beta1, dEOS, d2EOS
                if fkt == spf.sly_zeta:
                    para = [0.0, 8*math.pi, alfa, beta, spf.sly_dzeta, spf.sly_d2zeta] # rho (is only set to zero for initializing, in the next loop it will be set to its value and adapted everytime), kappa, alpha1, beta1, dEOS, d2EOS
                if fkt == spf.fps_zeta:
                    para = [0.0, 8*math.pi, alfa, beta, spf.fps_dzeta, spf.fps_d2zeta] # rho (is only set to zero for initializing, in the next loop it will be set to its value and adapted everytime), kappa, alpha1, beta1, dEOS, d2EOS
                
                # initialize the arrays I need for printing fR as a function of r (for each new integration they are initialized anew)
            
                while y0[0] > p_stop: # as long as the pressure is bigger than our final pressure, perform the integration
                    try:
                        #gsi = rf.rootfinding(fkt,np.log10(y0[0]/uconv[0]),0) # find the gsi associated with the pressure numerically
                        #rho = math.pow(10,gsi)*uconv[1] #calculate rho (I go into the TOV eq with it)
                        rho = qinterpol(y0[0])
                        para[0]= rho # save the rho into my para array (here the zero from before gets updated in first loop)
                        #f_R = 1.0+2.0*alfa*para[1]*(rho-3.0*y0[0]) # calculates fR (always before making the new step, i.e. the last values, i.e. r_final and fR_final will be missing in plot)
                        out = rk4.rk4_step(mod_tov_fRQ.modified_TOV, y0, r, h, para) # makes the RK4 step and saves the new pressure, mass, ideal stepsize and radious in the variable out
                        y0[0] = out[0] # saves new pressure
                        y0[1] = out[1] # saves new mass
                        #print(str(y0[0])+"   "+str(y0[1]))
                        r = r+h
                    except ValueError:
                        break
                    
                endtime = time.time() # calculates the time at the end of integration
            
                # print elapsed time to monitor integration speed
                print("Elapsed time in seconds")
                print(endtime-starttime)
                
                
                # print the final values
                print("Final Values")
                print(r/100000.0)
                print(y0[1]/(uconv[1]*1.99*math.pow(10,33)))
                
                # depending on which EOS I am using, append the results to different arrays
                if fkt == "HS_DD2.dat":
                    printimass_HSDD2.append(y0[1]/(uconv[1]*1.99*math.pow(10,33)))
                    printiradius_HSDD2.append(r/100000.0)
                    myfile.write(str(gsi_c)+" "+str(y0[1]/(uconv[1]*1.99*math.pow(10,33)))+" "+str(r/100000.0)+" HSDD2 "+str(alfa)+" "+str(beta)+"\n")
                if fkt == spf.sly_zeta:
                    printimass_sly.append(y0[1]/(uconv[1]*1.99*math.pow(10,33)))
                    printiradius_sly.append(r/100000.0)
                    myfile.write(str(gsi_c)+" "+str(y0[1]/(uconv[1]*1.99*math.pow(10,33)))+" "+str(r/100000.0)+" SLY "+str(alfa)+" "+str(beta)+"\n")
                if fkt == spf.fps_zeta:
                    printimass_fps.append(y0[1]/(uconv[1]*1.99*math.pow(10,33)))
                    printiradius_fps.append(r/100000.0)
                    myfile.write(str(gsi_c)+" "+str(y0[1]/(uconv[1]*1.99*math.pow(10,33)))+" "+str(r/100000.0)+" FPS "+str(alfa)+" "+str(beta)+"\n")
                
                printigsi.append(gsi_c) # append the gsi
                
                
                gsi_c = gsi_c+0.1 # increase gsi by a small value
                
                
            if fkt == "HS_DD2.dat":
                full_HSDD2.append([printimass_HSDD2, printiradius_HSDD2, alfa, beta])
            if fkt == spf.sly_zeta:
                full_sly.append([printimass_sly, printiradius_sly, alfa, beta])
            if fkt == spf.fps_zeta:
                full_fps.append([printimass_fps, printiradius_fps, alfa, beta])
    

myfile.close()

####################### GR CASE ############################

printimass_HSDD2_gr = []
printimass_sly_gr = []
printimass_fps_gr = []
printiradius_HSDD2_gr = []
printiradius_sly_gr = []
printiradius_fps_gr = []

gr_file = open("stardata_gr_georg.dat", "w")

tol = math.pow(10, -5) # sets the error tolerance between two rkf45 steps

for fkt in myfunctarray:
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

    gsi_c = 14.4 # set the value for the lowest gsi we are considering
    printigsi_gr = [] # initialize the array where we save all the gsi values we used in
    uconv = [p_conv, rho_conv] # put the conversion factors into an array
    gsi_final = 15.6
    # loop over all values of gsi until gsi = 16
    while gsi_c < gsi_final:
    
        m0 = 0 # initialize the mass (mass at the core of the neutron star is zero)
        #p_c = math.pow(10,fkt(gsi_c))*uconv[0] # calculate the central pressure
        p_c = pinterpol(math.pow(10, gsi_c)*uconv[1])
        t = 0 # set the integration parameter (will be the radius later) to zero
        #p_stop = p_c*math.pow(10,-12) # calculate the pressure when we stop the integration, we stop when the pressure is 12 orders of magnitude lower than the initial pressure
        p_stop = pressure[0]
        y0 = np.array([p_c, m0]) # array with initial values

        # print my starting values so that I can monitor progress
        print("Start Values")
        print(gsi_c)
        print(p_c)
        
        r = (t+0.0001)*h # start the integration process a little bit off the center in order to avoid the numerical singularity at r=0
        
        starttime = time.time() # calculate the starting time and save it
            
        # as long as the pressure is bigger than our final pressure, perform the integration    
        while y0[0] > p_stop:
            #gsi = rf.rootfinding(fkt,np.log10(y0[0]/uconv[0]),0) # find the gsi associated with the pressure numerically
            #rho = math.pow(10,gsi)*uconv[1] #calculate rho (I go into the TOV eq with it)
            rho = qinterpol(y0[0])
            #out = rkf.rk4_step(tov.TOV,y0,r,h,tol,rho) # makes the RKF45 step and saves the new pressure, mass, ideal stepsize and radious in the variable out
            out = rk4.rk4_step(tov.TOV, y0, r, h, rho)
            y0[0] = out[0]
            y0[1] = out[1]
            r = r+h
        
        endtime = time.time() # calculates the time at the end of integration
        
        # print elapsed time to monitor integration speed
        print("Elapsed time in seconds")
        print(endtime-starttime)
            

        # print the final values
        print("Final Values")
        print(r)
        print(y0[1])
        
        # depending on which EOS I am using, append the results to different arrays
        if fkt == "HS_DD2.dat":
            printimass_HSDD2_gr.append(y0[1]/(uconv[1]*1.99*math.pow(10,33)))
            printiradius_HSDD2_gr.append(r/100000.0)
            gr_file.write(str(gsi_c)+" "+str(y0[1]/(uconv[1]*1.99*math.pow(10,33)))+" "+str(r/100000.0)+" HSDD2 "+str(0)+" "+str(0)+"\n")
        if fkt == spf.sly_zeta:
            printimass_sly_gr.append(y0[1]/(uconv[1]*1.99*math.pow(10,33)))
            printiradius_sly_gr.append(r/100000.0)
            gr_file.write(str(gsi_c)+" "+str(y0[1]/(uconv[1]*1.99*math.pow(10,33)))+" "+str(r/100000.0)+" SLY "+str(0)+" "+str(0)+"\n")
        if fkt == spf.fps_zeta:
            printimass_fps_gr.append(y0[1]/(uconv[1]*1.99*math.pow(10,33)))
            printiradius_fps_gr.append(r/100000.0)
            gr_file.write(str(gsi_c)+" "+str(y0[1]/(uconv[1]*1.99*math.pow(10,33)))+" "+str(r/100000.0)+" FPS "+str(0)+" "+str(0)+"\n")
        
        printigsi_gr.append(gsi_c) # append the gsi
        
        
        gsi_c = gsi_c+0.1 # increase gsi by a small value

gr_file.close()
fulltime2 = time.time()
print("Total time needed in seconds")
print(fulltime2-fulltime1)
#pdb.set_trace()

#### Plots central pressure to mass ####
plt.figure(dpi=300)
for idx in range(len(full_HSDD2)):
    plt.plot(np.array(printigsi), np.array(full_HSDD2[idx][0]), label=r'$f(R,Q)$ with $\alpha =$'+str('{:.1e}'.format(full_HSDD2[idx][2]))+r' and $\beta=$'+str('{:.1e}'.format(full_HSDD2[idx][3])), linewidth = 0.5)
plt.plot(np.array(printigsi_gr), np.array(printimass_HSDD2_gr), label="GR", linewidth = 1, color="black")
plt.legend(prop={'size': 8}, title="HSDD2", loc="upper right")
plt.xlabel(r'$log_{10} \rho_c\ [\frac{g}{cm^3}]$')
plt.ylabel(r'M $[m_{\odot}]$')
plt.xlim(14,16.5)
plt.ylim(0,4)
plt.savefig("m_vs_rho_c_HSDD2_fRQ.png")
plt.clf()

#for idx in range(len(full_sly)):
#    plt.plot(np.array(printigsi), np.array(full_sly[idx][0]), label=r'$f(R,Q)$ with $\alpha =$'+str('{:.1e}'.format(full_sly[idx][2]))+r' and $\beta=$'+str('{:.1e}'.format(full_HSDD2[idx][3])), linewidth = 0.5)
#plt.plot(np.array(printigsi_gr), np.array(printimass_sly_gr), label="GR", linewidth = 1, color="black")
#plt.legend(prop={'size': 8}, title="SLy", loc="upper right")
#plt.xlabel(r'$log_{10} \rho_c\ [\frac{g}{cm^3}]$')
#plt.ylabel(r'M $[m_{\odot}]$')
#plt.xlim(14,16.5)
#plt.ylim(0,4)
#plt.savefig("m_vs_rho_c_sly_fRQ.png")
#plt.clf()

#for idx in range(len(full_fps)):
#    plt.plot(np.array(printigsi), np.array(full_fps[idx][0]), label=r'$f(R,Q)$ with $\alpha =$'+str('{:.1e}'.format(full_fps[idx][2]))+r' and $\beta=$'+str('{:.1e}'.format(full_HSDD2[idx][3])), linewidth = 0.5)
#plt.plot(np.array(printigsi_gr), np.array(printimass_fps_gr), label="GR", linewidth = 1, color="black")
#plt.legend(prop={'size': 8}, title="FPS", loc="upper right")
#plt.xlabel(r'$log_{10} \rho_c\ [\frac{g}{cm^3}]$')
#plt.ylabel(r'M $[m_{\odot}]$')
#plt.xlim(14,16.5)
#plt.ylim(0,4)
#plt.savefig("m_vs_rho_c_fps_fRQ.png")
#plt.clf()

#### End of Plotting central pressure to mass ####


#### Plots mass to radius ####
plt.figure(dpi=300)
for idx in range(len(full_HSDD2)):
    plt.plot(np.array(full_HSDD2[idx][1]), np.array(full_HSDD2[idx][0]), label=r'$f(R,Q)$ with $\alpha =$'+str('{:.1e}'.format(full_HSDD2[idx][2]))+r' and $\beta=$'+str('{:.1e}'.format(full_HSDD2[idx][3])), linewidth = 0.5)
plt.plot(np.array(printiradius_HSDD2_gr), np.array(printimass_HSDD2_gr), label="GR", linewidth = 1, color="black")
plt.legend(prop={'size': 8}, title="HSDD2", loc="upper right")
plt.xlabel(r'$R\ [km]$')
plt.ylabel(r'M $[m_{\odot}]$')
plt.xlim(8,20)
plt.ylim(0,4)
plt.savefig("m_vs_r_HSDD2_fRQ.png")
plt.clf()

#for idx in range(len(full_sly)):
#    plt.plot(np.array(full_sly[idx][1]), np.array(full_sly[idx][0]), label=r'$f(R,Q)$ with $\alpha =$'+str('{:.1e}'.format(full_sly[idx][2]))+r' and $\beta=$'+str('{:.1e}'.format(full_HSDD2[idx][3])), linewidth = 0.5)
#plt.plot(np.array(printiradius_sly_gr), np.array(printimass_sly_gr), label="GR", linewidth = 1, color="black")
#plt.legend(prop={'size': 8}, title="SLy", loc="upper right")
#plt.xlabel(r'$R\ [km]$')
#plt.ylabel(r'M $[m_{\odot}]$')
#plt.xlim(8,20)
#plt.ylim(0,4)
#plt.savefig("m_vs_r_sly_fRQ.png")
#plt.clf()

#for idx in range(len(full_fps)):
#    plt.plot(np.array(full_fps[idx][1]), np.array(full_fps[idx][0]), label=r'$f(R,Q)$ with $\alpha =$'+str('{:.1e}'.format(full_fps[idx][2]))+r' and $\beta=$'+str('{:.1e}'.format(full_HSDD2[idx][3])), linewidth = 0.5)
#plt.plot(np.array(printiradius_fps_gr), np.array(printimass_fps_gr), label="GR", linewidth = 1, color="black")
#plt.legend(prop={'size': 8}, title="FPS", loc="upper right")
#plt.xlabel(r'$R\ [km]$')
#plt.ylabel(r'M $[m_{\odot}]$')
#plt.xlim(8,20)
#plt.ylim(0,4)
#plt.savefig("m_vs_r_fps_fRQ.png")
#plt.clf()

#### End of plotting mass to radius ####

print("There's a starman waiting in the sky")
