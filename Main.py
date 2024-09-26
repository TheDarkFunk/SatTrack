#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 14:11:19 2019

     This file is a protoype for the final package main script
     It loads in a single TLE and station File and produces Ephemeris files for STK testing
     TO DO:
         Visibility
         Link Budget
         AOS/LOS Table
         User Input for TLE and Station
         Loop for all 32 satellites
@author: arvin
"""

from Fileio import banner,ReadStationFile,ReadNoradTLE,anykey, STKout
from datefun import ep2dat
from TrackingFun import TrackingTime, topo_to_range, PointingFile, max_range, LinkBudget,PointingFileSTK
from OMPython import ModelicaSystem
import matplotlib.pyplot as plt
#%% Initial User Input
banner()
anykey()
(t_start,t_stop,t_interval,t_step) = TrackingTime()
STNFIL = "station.dat"
print("##Reading Station File...")
Station = ReadStationFile(STNFIL)
print("##Reading TLE...")
f = open("TLE.txt","r")
Satellite=list()
sat_range = [[0]*1 for i in range(68)]
maximum_range = list()
p_sat_topo1=[[0]*1 for i in range(68)]
p_sat_topo2=[[0]*1 for i in range(68)]
p_sat_topo3=[[0]*1 for i in range(68)]
p_sat_ECF1=[[0]*1 for i in range(68)]
p_sat_ECF2=[[0]*1 for i in range(68)]
p_sat_ECF3=[[0]*1 for i in range(68)]
p_sat_ECI1=[[0]*1 for i in range(68)]
p_sat_ECI2=[[0]*1 for i in range(68)]
p_sat_ECI3=[[0]*1 for i in range(68)]
v_sat_topo1=[[0]*1 for i in range(68)]
v_sat_topo2=[[0]*1 for i in range(68)]
Rv_sat_topo3=[[0]*1 for i in range(68)]
v_sat_ECF1=[[0]*1 for i in range(68)]
v_sat_ECF2=[[0]*1 for i in range(68)]
v_sat_ECF3=[[0]*1 for i in range(68)]
v_sat_ECI1=[[0]*1 for i in range(68)]
v_sat_ECI2=[[0]*1 for i in range(68)]
v_sat_ECI3=[[0]*1 for i in range(68)]
elevation = [[0]*1 for i in range(68)]
azimuth = [[0]*1 for i in range(68)] 
eldot = [[0]*1 for i in range(68)]
azdot = [[0]*1 for i in range(68)] 
TableLine=[[0] * 4 for i in range(68)]
AOS = [[0]*1 for i in range(68)] 
LOS = [[0]*1 for i in range(68)]
gmst = [[0]*1 for i in range(68)]
time = [[0]*1 for i in range(68)]
power = [[0]*1 for i in range(68)]
R = [[0]*1 for i in range(68)]
doppler_shift = [[0]*1 for i in range(68)]
print("##Connecting to OpenModelica SatTrak Model")
#setup OM Server
#mod represents model
mod = ModelicaSystem("/home/arvin/Documents/Uni_4/ENG4350/4Tracking/OMdirectory/SatTrak.mo","SatTrak.SatTest", ["Modelica.Constants"])
for i in range(31):   
    line0 = f.readline()
    line1 = f.readline()
    line2 = f.readline()
    Satellite= ReadNoradTLE(line0, line1, line2)

#%% Simulation
    print("##### SATELLITE NO ",i+1," #####" )

    JD = float(line1[20:32])   + 2458484 - 2451545 - 0.5  # days from reference epoch to J2000
    simtime = t_start   + 2458484 - 2451545 - 0.5 #days from tracking day to J2000
    deltatime = simtime-JD #difference between TLE reference epoch and tracking time, in days
    newstarttime = deltatime*(86400) #converted to seconds
#input parameters...
    print("##Passing TLE parameters")
    mod.setParameters(**{'SatTestP1.M0':Satellite.meanan,'SatTestP1.N0':Satellite.meanmo,'SatTestP1.eccn':Satellite.eccn,'SatTestP1.Ndot2':Satellite.ndot,'SatTestP1.Nddot6':Satellite.nddot6,'SatTestP1.raan0':Satellite.raan,'SatTestP1.incl':Satellite.inc,'SatTestP1.argper0':Satellite.argper,'SatTestP1.tstart':newstarttime})
    print("##Passing Station parameters")
    mod.setParameters(**{'SatTestP2.stn_long':Station.stnlong,'SatTestP2.stn_lat':Station.stnlat,'SatTestP2.stn_elev':Station.stnalt})
    mod.setParameters(day=int(simtime),hour0=(simtime-int(simtime)+0.5))
    mod.getParameters()
#set simulation time
    mod.setSimulationOptions(startTime= 0, stopTime=t_interval, stepSize=t_step)
#Simulate Results
    print("##SIMULATING... \n...")
    mod.simulate()
    print("##SIMULATION COMPLETE")
    print("##RETRIEVING RESULTS")
#%% Process Results 
#
    print('...Time and GMST')
    (time[i],gmst[i]) = mod.getSolutions("time","GMST")
    print('...Look Angles')
    (azimuth[i],elevation[i]) = mod.getSolutions("azimuth","elevation")
    print('...Look Angle Rates of Change')
    (azdot[i],eldot[i]) = mod.getSolutions("azdot","eldot")
###(azdot,eldot) = mod.getSolutions("AzimuthRate","ElevationRate")
###get satellite trajectory for each coordinate system
###perifocal
#    print('...Peri Position')
#    (p_sat_pf1,p_sat_pf2,p_sat_pf3) = mod.getSolutions("p_sat_pf[1]","p_sat_pf[2]","p_sat_pf[3]")
#    print('...Peri Velocity')
#    (v_sat_pf1,v_sat_pf2,v_sat_pf3)= mod.getSolutions("v_sat_pf[1]","v_sat_pf[2]","v_sat_pf[3]")
###ECI
#    print('...ECI Position')
#    (p_sat_ECI1[i],p_sat_ECI2[i],p_sat_ECI3[i]) = mod.getSolutions("p_sat_ECI[1]","p_sat_ECI[2]","p_sat_ECI[3]")
#    print('...ECI Velocity')
#    (v_sat_ECI1[i],v_sat_ECI2[i],v_sat_ECI3[i]) = mod.getSolutions("v_sat_ECI[1]","v_sat_ECI[2]","v_sat_ECI[3]")
#
###ECF
#    print('...ECF Position')
#    (p_sat_ECF1[i],p_sat_ECF2[i],p_sat_ECF3[i]) = mod.getSolutions("p_sat_ECF[1]", "p_sat_ECF[2]", "p_sat_ECF[3]")
#    print('...ECF Velocity')
#    (v_sat_ECF1[i],v_sat_ECF2[i],v_sat_ECF3[i]) = mod.getSolutions("v_sat_ECF[1]", "v_sat_ECF[2]", "v_sat_ECF[3]")
###topo
    print('...Topo Position')
    (p_sat_topo1[i],p_sat_topo2[i],p_sat_topo3[i]) = mod.getSolutions("p_sat_topo[1]","p_sat_topo[2]","p_sat_topo[3]")
    print('...Range')
    R[i] = mod.getSolutions("R")
#    print('...Topo Velocity')
#    (v_sat_topo1[i],v_sat_topo2[i],v_sat_topo3[i]) = mod.getSolutions("v_sat_topo[1]","v_sat_topo[2]","v_sat_topo[3]")
#    print('...Visibility')
    (AOS[i],LOS[i]) = mod.getSolutions("AOS","LOS")

    print("##RESULTS RETRIEVED")
#plt.plot(time,p_sat_ECF1,time,p_sat_ECF2,time,p_sat_ECF3)
#%% Generate ephemeris files
#    refepoch = ep2dat(str(Satellite.refepoc+deltatime))
#    print("##WRITING EPHEMERIS FILES")
##    print("...'ephemPF.e'")
##    STKout("ephemPF", refepoch, time, "Custom PerifocalSystem Satellite/PRN_01_37753", p_sat_pf1, p_sat_pf2, p_sat_pf3, v_sat_pf1, v_sat_pf2, v_sat_pf3)
#    print("...'ephemECI.e'")
#    STKout("ephemECI"+str(i), refepoch, time[i], "Custom TEME CentralBody/Earth", p_sat_ECI1[i], p_sat_ECI2[i], p_sat_ECI3[i], v_sat_ECI1[i], v_sat_ECI2[i], v_sat_ECI3[i])
#    print("...'ephemECF.e'")
#    STKout("ephemECF"+str(i), refepoch, time[i], "Fixed", p_sat_ECF1[i], p_sat_ECF2[i], p_sat_ECF3[i], v_sat_ECF1[i], v_sat_ECF2[i], v_sat_ECF3[i])
#    print("...'ephemTopo.e'")
#    STKout("ephemTopo"+str(i), refepoch, time[i], "Custom TopocentricSystem Facility/Algonquin", p_sat_topo1[i], p_sat_topo2[i], p_sat_topo3[i], v_sat_topo1[i], v_sat_topo2[i], v_sat_topo3[i])
#    print("##EPHEMERIS FILES WRITTEN")
      
#%% Setting up Table
    
    print('##CALCULATING')
    print('...Ranges')
    #sat_range[i]=topo_to_range(p_sat_topo1[i],p_sat_topo2[i],p_sat_topo3[i])
    #maximum_range.append(max_range(sat_range[i],AOS[i],LOS[i],t_step))
    print('...Link Budgets')
    power[i] = LinkBudget(R[i])
    #power = -20#!! Placeholder
    TableLine[i][0]= Satellite.name
#converting AOS and LOS from simulation epoch seconds to the JulianDay format used by TLE, and then converting to UTC
#    TableLine[i][1] = ep2dat(str(t_start+ max(AOS[i])/86400 +19000))
#    TableLine[i][2] = ep2dat(str(t_stop+ max(LOS[i])/86400 +19000))
    TableLine[i][1] = 'None'
    TableLine[i][2] = 'None'
    
    if elevation[i][0] >= 9 and elevation[i][0] <= 89: #If starting inside FOR, set AOS to start time
        TableLine[i][1] = ep2dat(str(t_start+ max(AOS[i])/86400 +19000))        

    if max(LOS[i])>0.0: #if leaving FOR, set LOS
        TableLine[i][2] = ep2dat(str(t_start+ max(LOS[i])/86400 +19000))
        
    if max(AOS[i])>0.0: #if entering FOR, set AOS
        TableLine[i][1] = ep2dat(str(t_start+ max(AOS[i])/86400 +19000))
    if elevation[i][0]>=9 and elevation[i][0]<=89 and max(LOS[i])==0.0:
        TableLine[i][2]= ep2dat(str(t_stop+19000))
    if max(AOS[i])>0.0 and max(LOS[i])==0.0:
        TableLine[i][2] = ep2dat(str(t_stop+19000))
    TableLine[i][3] = min(power[i])

        
f.close()

#%% Write Visibility Table to User
print('AOS/LOS Table for ARO Facility. \n Computed from', ep2dat(str(t_start+19000)), 'to', ep2dat(str(t_stop+19000)), 'with timestep =', t_step, 'seconds \n')
print('{:2s} {:3s} {:28s} {:3s} {:30s} {:3s} {:30s} {:3s} {:15s}'.format(' #','|','Satellite Name','|','AOS','|','LOS','|','Min. Power [dBm]'))
print('')
for i in range(31):
    print('{:2d} {:3s} {:28s} {:3s} {:30s} {:3s} {:30s} {:3s} {:15f}'.format(i+1, '|', TableLine[i][0],'|', TableLine[i][1], '|', TableLine[i][2], '|', TableLine[i][3]))
anykey()
PointingFile(TableLine,AOS,LOS,azimuth,elevation,azdot,eldot,t_start,t_step)

#PointingFileSTK(TableLine,AOS,LOS,azimuth,elevation,azdot,eldot,t_start,t_step,time)