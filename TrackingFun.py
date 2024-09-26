#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 11:52:31 2019

Tracking Functions
 This python script contains some extra functions for the tracking project
      -Link Budget
      -Visibility
      -AOS/LOS Table generation
      -Output pointing file
@author: arvin
"""
import numpy as np
import datetime as dt
from Fileio import ERRMSG
from datefun import ep2dat2
def TrackingTime():
    today = dt.date.today()
    tt = today.timetuple()
    doy = tt.tm_yday
    
    print('You Will Be Prompted To Enter The Tracking Interval, in HH:MM:SS [UTC] \n     Example: 2:30 PM to 3 PM local time would be:  "Start: 18:30:00 Stop: 19:00:00"')
    start = input('Start Time: ')        
    bob = dt.timedelta(days=int(doy),hours=int(start[0:2]),minutes=int(start[3:5]),seconds=int(start[6:8]))
    t_start = bob.total_seconds()/86400 
    stop = input('Stop Time: ')
    jim = dt.timedelta(days=int(doy),hours=int(stop[0:2]),minutes=int(stop[3:5]),seconds=int(stop[6:8]))
    t_stop = jim.total_seconds()/86400
    t_interval = jim.total_seconds() - bob.total_seconds()
    if t_stop < t_start:
        ERRMSG('  Invalid Tracking Interval \n  Stop Time Must Be Greater Than Start Time')
        TrackingTime()
    else: 
        t_step = int(input('Timestep for simulation (seconds): ' ))
        print('## SIMULATING FOR ', dt.datetime.today().strftime('%Y-%m-%d'), 'from ', start, 'TO ', stop, ' [UTC]')
        return t_start,t_stop,t_interval,t_step
    
def LinkBudget(R):
    R = R*10**3
    C = list()
    Ls = list()
    f = 1575.42e6 #Frequency Band Center, Hz
    eta = 0.5 #Antenna Efficiency
    D = 46 #Antenna Diameter m
    c = 3e8 #speed of light, m/s
    pi = np.pi
    Gr = np.log10((eta/c**2)*(pi*f*D)**2)*10 #recieve gain [dBW]
    Gs = 4 #TBD @ARO [dBW]
    EIRP = (26.1+26.8)/2 #placeholder value I got from a paper [dBW]
    for i in range(len(R)):
            Ls.append(np.log10((c/(4*pi*R[i]*f))**2)*10) #free space path loss [dBW]
            C.append( Ls[i] + Gr + Gs + EIRP + 30)
    return C

def PointingFile(TableLine,AOS,LOS,Azimuth,Elevation,azdot,eldot,t_start,t_step):
    i = int(input('Please Select Which Satellite To Track: ' ))-1
    
    if TableLine[i][1] == str('None') and TableLine[i][2] ==str('None'):
     ERRMSG('  Invalid Satellite Selection \n  Please Select a Satellite with Valid Visibility')
     PointingFile(TableLine,AOS,LOS,Azimuth,Elevation,azdot,eldot,t_start,t_step)
    else:
        print('##Writing Pointing File for ',TableLine[i][0])
        print('...')
        
        if Elevation[i][0] >= 9 and Elevation[i][0] <= 89: #If starting inside FOR, set AOS to start time
            aos_index=0      
        if max(LOS[i])>0.0: #if leaving FOR, set LOS
            los_index=np.where(LOS[i]>0.0)[0][0]
        if max(AOS[i])>0.0: #if entering FOR, set AOS
            aos_index = np.where(max(AOS[i])>0.0)[0][0]
        if Elevation[i][0]>=9 and Elevation[i][0]<=89 and max(LOS[i])==0.0:
            los_index = len(Elevation[i])
        if max(AOS[i])>0.0 and max(LOS[i])==0.0:
            aos_index = np.where(AOS[i]>0)[0][0]
            los_index = len(Elevation[i])        
        f= open("Pointing.txt","w+") ## open a new text file
        f.write("{:19s} {:28s} \n".format("#Pointing File For ",TableLine[i][0]))
        #f.write("#---------------------------------------------------------------------- \n")
        #f.write('{:18s} {:25s} {:28s} \n'.format("# UTC Date/Time","Azimuth and AZ_Velocity", "Elevation and EL_Velocity"))
        #f.write('{:18s} {:25s} {:28s} \n'.format("#YYYY.DOY.HH:MM:SS", "AZd AZm AZ.s AZ.Vel", "ELd Elm EL.s EL.Vel"))
        for j in range(aos_index,los_index):
            AZd = int(Azimuth[i][j])
            AZm = int((Azimuth[i][j]-AZd)*60)
            AZs = ((Azimuth[i][j]-AZd)*60-AZm)*60
            AZv = azdot[i][j]
            
            ELd = int(Elevation[i][j])
            ELm = int((Elevation[i][j]-ELd)*60)
            ELs = ((Elevation[i][j]-ELd)*60-ELm)*60
            ELv = eldot[i][j]

            f.write('{:16s} {:3d} {:2d} {:4.1f} {:6.3f} {:2d} {:2d} {:4.1f} {:6.3f} \n'.format(ep2dat2(str(t_start+19000+j*t_step/86400)), AZd, AZm,AZs,AZv,ELd,ELm,ELs,ELv))    
        f.close()
        print('##COMPLETE')      
 
def PointingFileSTK(TableLine,AOS,LOS,Azimuth,Elevation,azdot,eldot,t_start,t_step,time):
    i = int(input('Please Select Which Satellite To Track: ' ))-1
    
    if TableLine[i][1] == str('None') and TableLine[i][2] ==str('None'):
     ERRMSG('  Invalid Satellite Selection \n  Please Select a Satellite with Valid Visibility')
     PointingFile(TableLine,AOS,LOS,Azimuth,Elevation,azdot,eldot,t_start,t_step)
    else:
        print('##Writing Pointing File for ',TableLine[i][0])
        print('...')
        
        if Elevation[i][0] >= 9 and Elevation[i][0] <= 89: #If starting inside FOR, set AOS to start time
            aos_index=0      
        if max(LOS[i])>0.0: #if leaving FOR, set LOS
            los_index = np.where(max(LOS[i]))[0][0]
        if max(AOS[i])>0.0: #if entering FOR, set AOS
            aos_index = np.where(max(AOS[i]))[0][0]
        if Elevation[i][0]>=9 and Elevation[i][0]<=89 and max(LOS[i])==0.0:
            los_index = len(Elevation[i])
        f= open("PointingSTK.sp","w+") ## open a new text file
        f.write("stk.v.4.3 \n")
        f.write("Begin Attitude \n")
        f.write('NumberofAttitudePoints ')
        f.write(str(los_index-aos_index))
        f.write('\n')
        f.write('Sequence 323 \n')
        f.write('AttitudeTimeEzElAngles \n')
        for j in range(aos_index,los_index):
            f.write('{:8.3f} {} {} \n'.format(time[i][j],Azimuth[i][j],Elevation[i][j]))    
        f.write('End Attitude')
        f.close()
        print('##COMPLETE')      
              
def topo_to_range(p_sat_topo1,p_sat_topo2,p_sat_topo3):
    
    sat_range=list()
    for j in range(len(p_sat_topo1)):
        sat_range.append(np.sqrt(p_sat_topo1[j]**2+p_sat_topo2[j]**2+p_sat_topo3[j]**2))
    return sat_range

def max_range(sat_range,AOS,LOS,t_step):
    aos_index = int(round(max(AOS)/t_step))
    los_index = int(round(max(LOS)/t_step))
    maximum_range=list()
    if max(LOS)==0.0 or max(AOS)==0.0:
            maximum_range.append(max(sat_range))
    else:
            maximum_range.append(max(sat_range[aos_index:los_index]))
    return maximum_range