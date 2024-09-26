#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 13:39:44 2019
Updated March 21st 2019

FileIO.py
all the I/O needed for front end user information
uses the functions:
       Banner()
       anykey()
       errmsg()
       ReadStationFile()
       ReadNoradTLE()
       STKout()
@author: arvin
"""
#importing
import datetime
from colorama import Fore, Style
from pprint import pprint
from collections import namedtuple
#%%
#Defining Functions
def banner():
    """
Created on Thu Feb 28 13:13:17 2019
Banner()
Prints the ENG4350 team info, program name, version and welcome message
@author: arvin
"""
    vernum = 'Version 3.0'
    print('-------SatTrak------- \n Arvin T. & Tiffany J. \n', datetime.datetime.now(), '\n', vernum, ' \n -------WELCOME-------')
#%%
def anykey():
    input('Press Enter to Continue: ')
#%%    
def ERRMSG(STRING):
    """
Created on Thu Feb 28 13:33:16 2019
ERRMSG(STRING)
"beeps" the terminal, displays an error message STRING and exits
@author: arvin
"""
    
    print(Fore.RED + 'ERROR')
    print(STRING)
    print(Style.RESET_ALL)
#%%
def ReadStationFile(STNFIL):
    """
Created on Fri Mar  1 11:43:12 2019
ReadStationFile(STNFIL)
Returns a Python “named tuple” with the parameters defined below. This function reads
the station parameters file and returns the values in degrees and meters for geographical
coordinates, degrees for az/el limits and deg/min for maximum az/el speeds.
Format for STNFIL.dat file:
Station name (any string)
Latitude (float, degrees)
Longitude (float, degrees. Negative west of Greenwich)
Altitude (meters)
Local time zone shift with daylight saving hour and proper sign, Eastern time -5 or -4.
@author: arvin
"""
    #open file
    f = open(STNFIL,"r");
    mylist = f.read().splitlines()#read all parameters as a lit
    #create a named tuple for station info
    classname = 'Station'
    fields = 'name stnlat stnlong stnalt utc_offset az_el_nlim az_el_lim st_az_speed_max st_el_speed_max'
    Station = namedtuple(classname, fields)
    #populate our namedtuple with station.dat information
    Station.name = mylist[0];
    Station.stnlat = float(mylist[1]);
    Station.stnlong = float(mylist[2]);
    Station.stnalt= float(mylist[3]);
    Station.utc_offset = float(mylist[4]);
    Station.az_el_nlim = int(mylist[5]);#number of azimuth restrictions
    #creating a new named tuple to store the min and max elevation at restricted azimuth
    Station.az_el_lim = namedtuple('az_el_lim', 'az elmin elmax');
    #populate this tuple with values, accoridng to number of restrictions
    
    Station.az_el_lim.az = float(mylist[6][0:3])
    Station.az_el_lim.elmin = float(mylist[6][4:8])
    Station.az_el_lim.elmax = float(mylist[6][10:14])
    
    Station.st_az_speed_max = float(mylist[7]);
    Station.st_el_speed_max = float(mylist[8]);
    f.close()
    #out = (Station.name,Station.stnlat,Station.stnlong,Station.stnalt,Station.utc_offset,Station.az_el_nlim, Station.az_el_lim.az,Station.az_el_lim.elmax,Station.az_el_lim.elmin, Station.st_az_speed_max,Station.st_el_speed_max);
    return Station #pprint(dict(vars(Station)))
   
#%%    
def ReadNoradTLE(line0,line1,line2):
    """
Created on Sat Mar  2 16:31:42 2019
Satellite[i] = ReadNoradTLE(line0, line1, line2)
Satellite[i] is a Python list of named tuples, each tuple has elements
Satellite[i].name
Satellite[i].refepoch,
Satellite[i].incl,
Satellite[i].raan,
Satellite[i].eccn,
Satellite[i].argper,
Satellite[i].meanan,
Satellite[i].meanmo,
Satellite[i].ndot,
Satellite[i].nddot6,
Satellite[i].bstar,
Satellite[i].orbitnum
This function reads 3 lines (the last two are in TLE format) and determines orbital
variables for the satellite.
@author: arvin
"""
    classname = 'Satellite'
    fields = 'name refepoc inc raan eccn argper meanan meanmo ndot nddot6 bstar orbitnum'
    Satellite= namedtuple(classname, fields)
    #getting orbital paremeters based on their index within line strings
    Satellite.name = line0[0:22];
    Satellite.refepoc = float(line1[18:32])
    Satellite.inc = float(line2[8:17])
    Satellite.raan = float(line2[17:26])
    Satellite.eccn = float('0.'+line2[26:33])#append 0. in front of entry for eccentricity
    Satellite.argper = float(line2[34:43])
    Satellite.meanan = float(line2[43:52])
    Satellite.meanmo = float(line2[52:64])
    Satellite.ndot = float(line1[33:44])
    Satellite.nddot6 = float(line1[44:50] + 'e' + line1[51:52])
    Satellite.bstar = float(line1[53:59] + 'e' + line1[60:62])
    Satellite.orbitnum = int(line2[63:69])
    #out = (Satellite.name,Satellite.refepoc,Satellite.inc,Satellite.raan,Satellite.eccn,Satellite.argper,Satellite.meanan,Satellite.meanmo,Satellite.ndot,Satellite.nddot6,Satellite.bstar,Satellite.orbitnum)
    return Satellite #pprint(dict(vars(Satellite)))

#%%
def STKout(outfile, StartEpoch, time, Coord, px, py, pz, vx, vy, vz):
    
   f= open(outfile +".e","w+") ## open a new text file
   f.write("stk.v.4.3 \n")
   f.write("BEGIN Ephemeris \n")
   f.write("NumberOfEphemerisPoints %d \n" % len(px))
   # read from OMPython scenepoch = input('What is the scenario Epoch?\n') 
   f.write("ScenarioEpoch " + StartEpoch + "\n")
   f.write("InterpolationMethod Lagrange \n")
   f.write("InterpolationOrder 7 \n")
   f.write("CentralBody Earth \n" )
   f.write("CoordinateSystem " + Coord + "\n")
   f.write("DistanceUnit Kilometers \n") 
   f.write("EphemerisTimePosVel \n \n")
   
   ## Insert data printing here
   for j in range(0,len(px)):
       f.write("%e %e %e %e %e %e %e \n" %  (time[j], px[j], py[j], pz[j], vx[j], vy[j], vz[j]))
   
   f.write("END Ephemeris")
   f.close()