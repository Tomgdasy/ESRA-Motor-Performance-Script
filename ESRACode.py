# Title: ESRA Motor Performance Calculator
#
# Author: Tom Gerendasy
#
# Last Updated: November, 25th 2020
#
# Description: This performance program will return delivered ISP, characteristic 
#              velocity, total impulse, coefficicent of thrust, and nozzle exit 
#              pressure. It will also return theoretical ISP, coefficient of 
#              thrust, and characteristic velocity
#
# Reminder: Make sure filenames are formatted properly!
from tkinter import filedialog
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from pygasflow.solvers import isentropic_solver
import math

#Importing CSV and compiling data
def data_collect():
    # Select one or more files
    # Files must have the format "~x.xxxxx~x.xxxxx~x.xxxxx.csv"
    # Where each place holder represents: "~Nozzle throat diameter [in.]~Nozzle exit diameter[in.]~mass propellant [lbs.].csv"
    import_file_path = filedialog.askopenfilenames()                           
    num = len(import_file_path) 
   
    # Variable allocation 
    L = np.ones(num).astype(int)                                                         
    Nd = np.ones(num).astype(float)
    Nde = np.ones(num).astype(float)
    mp = np.ones(num).astype(float)
    pMax = np.ones(num).astype(float)
    tMax = np.ones(num).astype(float)
    pIndex = np.ones(num).astype(int)
    tIndex = np.ones(num).astype(int)
    indexIS = np.ones(num).astype(int)
    indexIE = np.ones(num).astype(int)
    burnTime = np.ones(num).astype(float)   
    Bt = np.ones(num).astype(float)                              
    Bs = np.ones(num).astype(int)                                                      
    Be = np.ones(num).astype(int)
    timescale = 10**(-6)
    # Variable allocation end
   
    # Loop used to acquire motor information: Nozzle geometries, propellant mass
    i = 0  
    for y in import_file_path:
        L[i] = len(pd.read_csv(y))                                                 
        geom = []
        i2 = 0
        for z in y:
            if z == '~':
                geom.append(float(y[i2 + 1:i2 + 8]))
            i2 += 1
        Nd[i] = geom[0]
        Nde[i] = geom[1]
        mp[i] = geom[2]
        i += 1   
                                                            
    # data matrix variable allocation    
    data = np.zeros((max(L), num*3))                                                                                                                                                                
                                   
    # For loop used to read and record csvs into data matrix using pandas package    
    i = 0                                                                      
    for y in list(range(0,num*3,3)):                                         
         data_temp= pd.read_csv(import_file_path[i]).to_numpy()    
         # Recording in data matrix and eliminating potential bad data at beginning and end   
         data[0:L[i]-100,y:y+3] = data_temp[50:L[i]-50,0:3]
         data[:,y] = data[:,y] *timescale
         # Calibration for pressure transducer
         data[:,y+1] = data[:,y+1] * 0.064271748 - 750
         # Calibration for load cell
         data[:,y+2] = data[:,y+2] * 0.1628 - 262.98  
         
         [pMax[i], pIndex[i]] =[max(data[:,y+1]), np.argmax(data[:,y+1])]  
         tIndex[i] = np.argmax(data[:,y+2])
         # "Zeroing the scale" and converting to [lbs.]
         data[:,y+2] = (data[:,y+2] - np.mean(data[tIndex[i]-1000:tIndex[i]-500,y+2])) *  0.224809 
         tMax[i] = max(data[:,y+2]) 
       
         # Loop used to locate indices of burn start and end                                                                 
         for u in list(range(pIndex[i], 0, -1)):   
             # If pressure is less than 15 Psi, mark the index and break. Index Start.                                  
             if data[u,y+1] < 15:                                                 
                 indexIS[i] = u                                                     
                 break
         for u in list(range(pIndex[i], len(data))):
             # If pressure is less than 15 Psi, mark the index and break. Index End.                                 
             if data[u, y+1] < 15:                                                 
                 indexIE[i] = u                                             
                 break
         burnTime[i] = data[indexIE[i], y] - data[indexIS[i], y]
         
         # Loop used to locate indices of steady state burn start and end 
         for u in list(range(pIndex[i], 0, -1)):   
         # If pressure is less than 60 percent maximum pressure, index start.                                 
            if data[u,y + 1] < 0.6*pMax[i]:                                                 
             Bs[i] = u                                                     
             break
         for u in list(range(pIndex[i], len(data))):
         # # If pressure is less than 90 percent maximum pressure, index end.                                     
            if data[u, y + 1] < 0.9*pMax[i]:                                                 
             Be[i] = u                                             
             break
         Bt[i]  = data[Be[i],y] - data[Bs[i],y] 
         i += 1  


         
    return [data,num,pMax,pIndex,tMax,tIndex,indexIE,indexIS,burnTime,import_file_path,Nd,Nde,mp,Bs,Be,Bt] 
 
                                                 
[data,num,pMax,pIndex,tMax,tIndex,indexIE,indexIS,burnTime,import_file_path,Nd,Nde,mp,Bs,Be,Bt] = data_collect()

def performance():
    # Variable allocation
    impulse = np.ones(num).astype(float)
    impulse_ss = np.ones(num).astype(float)
    Isp = np.ones(num).astype(float)
    Isp_ss = np.ones(num).astype(float)
    Cstar = np.ones(num).astype(float)
    Cf = np.ones(num).astype(float)
    
    # Loop used to calculate impulse, ISP, characteristic velocity, and coefficient of thrust                                                             
    i = 0
    for y in list(range(0, num*3, 3)):     
        # Impulse [N-sec]. Calculating area under the thrust curve. Conversion to N-sec                                          
        impulse[i] = np.trapz(data[indexIS[i]:indexIE[i], y+2],                     
                              data[indexIS[i]:indexIE[i], y]) / 0.224809 

        # Specific Impulse [sec]. Dividing the impulse by the weight of the propellant. Propellant weight converted to N. 
        Isp[i] = impulse[i]/(mp[i]*4.448)

        # Characteristic Velocity [ft/s].Calculating area under the pressure curve. Propellant weight converted to lbm                                         
        Cstar[i]=(np.trapz(data[indexIS[i]:indexIE[i], y+1],                        
                          data[indexIS[i]:indexIE[i],
                                y])) * (np.pi/4*Nd[i]**2) / (mp[i] / 32.2)  
        # Coefficient of Thrust calculation. Maximum pressure and thrust                    
        Cf[i] = np.mean(data[Bs[i]:Be[i],y+2]/(data[Bs[i]:Be[i],y+1]*np.pi/4*Nd[i]**2))                                               
        i += 1
    return [impulse, Isp, impulse_ss, Isp_ss, Cf, Cstar]

[impulse, Isp, impulse_ss, Isp_ss, Cf, Cstar] = performance()

def plot():
    i = 0
    for y in list(range(0, num*3, 3)):
        plt.figure()
        # Plot Data Set                                                
        plt.plot(data[indexIS[i]-50:indexIE[i]+50, y], data[indexIS[i]-50:indexIE[i]+50, y + 1], color = 'k')
        # Plot peak pressure
        plt.plot(data[pIndex[i], y], data[pIndex[i], y+1], marker='*', color ='b', label = 'Max Pressure [' + str(round(pMax[i], 2))+' Psi]')
        # Plot steady state start time
        plt.plot(data[indexIS[i], y], data[indexIS[i], y+1], marker = '*', color = 'r',label = 'Burn Start')
        # Plot steady state end time
        plt.plot(data[indexIE[i], y], data[indexIE[i], y+1], marker = '*', color = 'g', label = 'Burn End')
        plt.grid(color='k', linestyle='-', linewidth=.12)
        plt.xlabel('Time [s]') 
        plt.ylabel('Combustion Pressure [psi]') 
        plt.title('Subscale Static Test, Nozzle Throat Diameter: ' + str(Nd[i]) + ' In.') 
        plt.legend(bbox_to_anchor=(.8,.58))
        
        
        plt.figure()
        # Plot Data Set                                                
        plt.plot(data[indexIS[i]-50:indexIE[i]+50, y], data[indexIS[i]-50:indexIE[i]+50, y + 2], color = 'dimgrey')
        # Plot peak pressure
        plt.plot(data[pIndex[i], y], data[pIndex[i], y+2], marker='*', color = 'b', label = 'Max Thrust [' + str(round(tMax[i], 2))+' lbs.]')
       
        plt.plot(data[indexIS[i], y], data[indexIS[i], y+2], marker = '*', color = 'r',label = 'Start and end times')
        plt.plot(data[indexIE[i], y], data[indexIE[i], y+2], marker = '*', color = 'r')
       
        plt.grid(color='k', linestyle='-', linewidth=.12)
        plt.xlabel('Time [s]') 
        plt.ylabel('Thrust [lbs.]') 
        plt.title('Subscale Static Test, Nozzle Throat Diameter: ' + str(Nd[i]) + ' In.') 
        plt.legend(bbox_to_anchor=(.8,.58))
        i += 1
    
    from tabulate import tabulate
    datatable= np.zeros((num, 6))    
    header = ['Motor', 'Max Pressure [Psi]', 'Total Impulse [N-sec]', 'Delivered ISP [sec]', 'C* [ft/sec]', 'Cf']
    for y in range(num):
        aa = [y+1,round(pMax[y], 3), round(impulse[y], 3),round(Isp[y], 3), round(Cstar[y], 3), round(Cf[y], 3)]
        datatable[y, 0] = aa[0]
        datatable[y, 1] = aa[1]
        datatable[y, 2] = aa[2]
        datatable[y, 3] = aa[3]
        datatable[y, 4] = aa[4]
        datatable[y, 5] = aa[5]
    print(tabulate(datatable, headers = header))

plot()

def nozzlePerformance():
    
    ppt3 = np.ones(num)
    
    for i in list(range(num)):
        Exp = Nde[i]**2*np.pi/4/(Nd[i]**2*np.pi/4)
        ppt3[i] =list(isentropic_solver('crit_area_super', Exp, 1.2016))[1]
        

    pNozzle = np.ones([len(data),num])
    pNozzleMax = np.ones(num)
    pNozzleIndex = np.ones(num)
    i=0
    # Loop used to calculate nozzle exit pressure using supersonic isentropic flow tables
    for y in list(range(0,num*3,3)):
        pNozzle[0:len(data[indexIS[i]-50:indexIE[i]+25,y+1]),i] = ppt3[i] * data[indexIS[i]-50:indexIE[i]+25,y+1]
        [pNozzleMax[i], pNozzleIndex[i]] =[max(pNozzle[:,i]), np.argmax(pNozzle[:,i])]
        
    
        plt.figure()
        tttarray = np.linspace(0,burnTime[i],len(pNozzle[0:len(data[indexIS[i]-50:indexIE[i]+25,y+1]),i]))
        plt.plot(tttarray,pNozzle[0:len(data[indexIS[i]-50:indexIE[i]+25,y+1]),i],'-', color = 'dimgray',
                  label = 'ISP = ' +str(round(Isp[i],2)) + ' Cf = ' + str(round(Cf[i],3)) + ' C* = ' + str(round(Cstar[i],2)))
        
        plt.plot(tttarray[pNozzleIndex[i].astype(int)],pNozzleMax[i],'*k')
        plt.plot(np.linspace(-1,3,10),14.7*np.ones(10),color = 'k', linestyle='--')
        
        plt.grid(color='k', linestyle='-', linewidth=.12)
        plt.xlabel('Time [s]') 
        plt.ylabel('Nozzle Exit Pressure [psi]') 
        plt.title('Nozzle Exit Pressure Vs. Time') 
        plt.legend(bbox_to_anchor=(1.4,.58))
        plt.xlim(0,burnTime[i])
        plt.ylim(-2,18)
        i=i+1
    return [pNozzleMax, pNozzleIndex]
[Pe, pNozzleIndex] = nozzlePerformance()


def theoreticals():
    #Calculating Correction factors

    nozzle_div = 0.5*(math.cos(15*np.pi/180)+1)
    p_corr = np.ones(num)
    comb_corr = np.ones(num)
    Isp_ideal = np.ones(num)
    Isp_theo = np.ones(num)
    Cstar_theo = np.ones(num)
    Cf_theo = np.ones(num)
    # Universal gas constant
    R = 8134
    # Combustion temp for LS Propellant
    T0 = 2538.31
    # 2 phase flow specific heat ratio for LS propellant
    k = 1.2016
    # Specific heat of mixture of LS propellant
    Kmix = 1.23528
    # Molecular weight of LS propellant
    Mw = 23.0847
    # Ambient Pressure
    Pa = 14.7
    i = 0
    # Loop used to calculate the pressure correction factor
    for y in range(0,num*3,3):
        p1 = np.trapz(data[indexIS[i]:Bs[i],y+1], data[indexIS[i]:Bs[i],y])
        p2 = np.trapz(data[indexIE[i]:Be[i],y+1], data[indexIE[i]:Be[i],y])
        p3 = np.trapz(data[indexIS[i]:indexIE[i],y+1], data[indexIS[i]:indexIE[i],y])
        p_corr[i] =1 - (p1+ p2)/p3
        i += 1
    for i in range(num):
        Isp_ideal[i]=1/9.806*(2*k/(k-1)*(R*T0/Mw)*(1-(Pa/pMax[i])**((k-1)/k)))**.5
        Cstar_theo = (R*T0/(Mw*Kmix)*((Kmix+1)/2)**((Kmix+1)/(Kmix-1)))**.5 * 3.28084
        Cf_theo[i] = k * (2/(k-1)*(2/(k+1))**((k+1)/(k-1))*(1-(Pa/pMax[i])**((k-1)/k)))**.5 
        
        comb_corr[i] = Cstar[i]/Cstar_theo
        Isp_theo[i] = comb_corr[i] * p_corr[i] * nozzle_div * 0.95 * Isp_ideal[i]
    return [Isp_theo,Cstar_theo,Cf_theo,k,R,Mw,Pa,Kmix,T0,Isp_ideal,p_corr,comb_corr] 
[Isp_theo,Cstar_theo,Cf_theo,k,R,Mw,Pa,Kmix,T0,Isp_ideal,p_corr,comb_corr] = theoreticals()
