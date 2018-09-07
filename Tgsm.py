import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import pandas as pd
import os
from matplotlib import cm
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
from power_to_temperature import Radio_source_trans
Cmap = cm.jet
Cmap.set_under("w")

def trajectory(time,lon = -118.3011,lat = 28.9733):
    """
    Calculates trayectories for the antenna in a given time.
    Returns an array of l and b galactic coordinates values in degrees, for a given location.
    
    It is asummed that the antenna is looking derectly to the zenit, 
    so in this case the right ascention (RA) is equal to the local sideral time of the location 
    and the declination (DEC) is equal to the latitude.
    
    Parameters:
    time: local time of observation, can be an array or a single time. Prefered format is
          'yyyy-mm-dd hh:mm:ss'.
          
           WARNING: Make sure that time is in UTC.
    
    Optional parameters:
    lon: Default longitude is given for Isla de Guadalupe at -118.3011 degrees
    lat: Default latitude is given for Isla de Guadalupe at 28.9733 degrees
    """
    t = Time(time, location =(lon,lat))
    RA = np.array(t.sidereal_time('mean').degree)
    DEC = lat*np.ones(len(RA))
    coords = SkyCoord(ra=RA*u.degree,dec=DEC*u.degree)
    l = coords.galactic.l.degree
    b = coords.galactic.b.degree
    return l,b

def pattern(time,Freq,PATH="antenna_beam/",lon = -118.3011,lat = 28.9733):
    """
    Calculates the beam pattern for the antenna at a given time.
    Returns an array of l and b galactic coordinates values in degrees of the antenna pattern, 
    for a given location and day.
    
    It also returns an array of Temperature lecture of the antenna for a given coordinate (l,b)
    
    Parameters:
    time: Time of observation, can be an array or a single time. Prefered format is
          'yyyy-mm-dd hh:mm:ss'. 
          
           WARNING: Make sure that time is in UTC.
          
    Freq: Frequency of the antenna beam in MHz.
    
    Optional parameters:
    PATH: Folder where the pattern is stored, note that the files within this folder
          must be named with its frequency, for example 70MHz.hdf5.
         
         Default is antenna_beam.
         
    lon: Default longitude is given for Isla de Guadalupe at -118.3011 degrees
    lat: Default latitud is given for Isla de Guadalupe at 28.9733 degrees
    """
    Data = pd.read_hdf(PATH+"0%dMHz.hdf5"%Freq) #Change path if beam pattern is changed
    t = Time(time, location =(lon,lat))
    LST = t.sidereal_time('mean').degree
    theta,phi = np.radians(Data.values[:,0]),np.radians(Data.values[:,1])
    dB = Data.values[:,2]
    X,Y,Z=np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta) 
    Rxy = np.sqrt(X**2.+Y**2.)
    colat,Az = np.arctan2(Rxy,Z),np.arctan2(Y,X)
    Alt = 0.5*np.pi-colat
    lat = np.radians(lat)
    sinDEC = np.sin(Alt)*np.sin(lat)-np.cos(Alt)*np.cos(Az)*np.cos(lat)
    DEC = np.arcsin(sinDEC)
    sinH = np.sin(Az)*np.cos(Alt)/np.cos(DEC)
    cosH = (np.cos(Az)*np.cos(Alt)*np.sin(lat)+np.sin(Alt)*np.cos(lat))/np.cos(DEC)
    H = np.arctan2(sinH,cosH)
    DEC = np.degrees(DEC)
    RA = LST - np.degrees(H)
    coords = SkyCoord(ra=RA*u.degree,dec=DEC*u.degree)
    l = coords.galactic.l.degree
    b = coords.galactic.b.degree
    Temp = Radio_source_trans(dB,Freq,1e6)
    return l,b,Temp

def convolve(time,Freq,PATH="antenna_beam/"):
    """
    Convolves the antenna beam pattern with the gsm map of the galaxy for a given frequenacy.
    Returns the convolved temperature of the gsm.
    
    Parameters:
    time: local time of observation, can be an array or a single time. Prefered format is
        'yyyy-mm-dd hh:mm:ss'.
        
        
         WARNING: Make sure that time is in UTC.
    Optional parameters:
    PATH: Folder where the pattern is stored, note that the files within this folder
          must be named with its frequency, for example 70MHz.hdf5.
         
         Default is antenna_beam.         
         
    Freq: Frequency desired in MHz.
    
    Optional parameters:
    PATH: Folder where the pattern is stored, note that the files within this folder
          must be named with its frequency, for example 70MHz.hdf5.
         
         Default is antenna_beam.
    """
    nside = 32
    Data = pd.read_hdf("gsm_maps/gsm_%dMHz.hdf5"%Freq)
    bmap_gal = Data.values[:,0]
    bmap_gal2 = hp.ud_grade(bmap_gal,nside)
    l,b,Temp = pattern(time,Freq,PATH)
    pix = hp.ang2pix(nside,l, b, lonlat=True)
    bmap_pat = np.zeros(hp.nside2npix(nside))
    bmap_pat[pix] = Temp
    T_gsm = sum(bmap_gal2*bmap_pat)/sum(bmap_pat)
    return T_gsm

def T_gsm(time,freqs=(50,90),bins=20,days=1,PATH="antenna_beam/"):
    """
    Provides a table of the convolved temperature of the GSM map with the 
    Antenna beam pattern for a full day of observation, in a range of frequencies.
    It saves a file named Tgsm.hdf5 with the values obtained  in the calibration folder.
    
    Parameters:
    time: Initial date and hour of observation, prefered format is 'yyyy-mm-dd hh:mm:ss'. 
        
            WARNING: Make sure that time is in UTC.
    
    Optional parameters:
    PATH: Folder where the pattern is stored, note that the files within this folder
          must be named with its frequency, for example 70MHz.hdf5.
         
         Default is antenna_beam.
    freqs: Range of frequencies, must be a tuple with initial frequency and final frequency.
           Default is 50-90
    
    bins: Observation interval in minutes, default is 20 minutes.
    days: Days of observation, default is 1 day.
    """
    if not os.path.exists('calibration'):
        os.makedirs('calibration')
    Freqs = np.arange(freqs[0],freqs[1]+1)
    t0 = Time(time)
    dt = bins*u.min # Modify unit of time interval if needed
    DT = dt.to(u.hour)
    times = t0 + DT*np.arange(0,days*24/DT.value)
    data = np.zeros([len(Freqs),len(times)])
    i,j = 0,0
    for f in Freqs:
        for k in range(len(times)):
            data[i,j] = convolve(times[k],f,PATH)
            j+=1
        i+=1
        j=0
    df = pd.DataFrame(data,index = Freqs,columns = times.value)
    df.to_hdf('calibration/Tgsm.hdf5','df')
    return df

def check_LST(time,lon = -118.3011,lat = 28.9733):
    """
    Checks the Local Sidereal time for a given time in a given location.
    
    Parameters:
    time: Time of observation, can be an array or a single time. Prefered format is
          'yyyy-mm-dd hh:mm:ss'. 
          
           WARNING: Make sure that time is in UTC.
          
    Optional parameters:
    lon: Default longitude is given for Isla de Guadalupe at -118.3011 degrees
    lat: Default latitude is given for Isla de Guadalupe at 28.9733 degrees
    """
    t = Time(time, location =(lon,lat))
    LST = t.sidereal_time('mean')
    print 'LST time:',LST
    return LST