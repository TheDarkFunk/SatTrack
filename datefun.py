# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 13:15:48 2019
@authors: Arvin Tangestanian and Tiffany Joseph
"""

## Define function doy(). Doy gives the number of days into the year of the current day. 
# This function enables the user to enter a year, month and date numerically, then outputs the Gregorian Calendar date.
def doy(YR,MO,D):

    month_days = [31,28,31,30,31,30,31,31,30,31,30,31]
    MO_p = MO-1
    if YR%4==0:
        month_days[1] = 29
    else:  month_days[1] = 28
    numberDays = sum(month_days[0:MO_p])+D
    return numberDays


## Define function dayFraction(). dayFraction() gives the fraction of the day that has elapsed
    #This function takes the time of day and calculates the fraction of the day at that time
def frcofd(HR,MI,SE):

    dayFraction = ((SE/60 + MI)/60 + HR)/24
    return dayFraction

def ep2dat(JulianDate): #input as a stirng
    import datetime as dt
    import numpy as np
    # Extract Reference Epoch
    refepochFloat  = float(JulianDate)
    # Extract the reference epoch year, month and day
    refepochdt = dt.datetime.strptime(JulianDate[:5],'%y%j')
    # Extract the reference epoch hour, minute and second
    dfrac = np.modf(refepochFloat)[0]
    [rem, hr] = np.modf(dfrac*24)
    [rem2, mins] = np.modf(60*rem)
    [rem3, secs] = np.modf(60*rem2)
    mics = np.modf(rem3*10**6)[1]
    refepochdt = refepochdt.replace(hour=np.int(hr), minute=np.int(mins), 
                                    second=np.int(secs), 
                                    microsecond=np.int(mics))
    return refepochdt.strftime('%d %b %Y %H:%M:%S.%f')

def ep2dat2(JulianDate): #input as a stirng
    import datetime as dt
    import numpy as np
    # Extract Reference Epoch
    refepochFloat  = float(JulianDate)
    # Extract the reference epoch year, month and day
    refepochdt = dt.datetime.strptime(JulianDate[:5],'%y%j')
    # Extract the reference epoch hour, minute and second
    dfrac = np.modf(refepochFloat)[0]
    [rem, hr] = np.modf(dfrac*24)
    [rem2, mins] = np.modf(60*rem)
    [rem3, secs] = np.modf(60*rem2)
    mics = np.modf(rem3*10**6)[1]
    refepochdt = refepochdt.replace(hour=np.int(hr), minute=np.int(mins), 
                                    second=np.int(secs), 
                                    microsecond=np.int(mics))
    #return refepochdt.strftime('%d %b %Y %H:%M:%S.%f')
    return refepochdt.strftime('%Y.%j.%H:%M:%S')