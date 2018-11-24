import pandas as pd
import numpy as np
import os
import glob
from power_to_temperature import *
from test import *

Freqs = np.linspace(1e-22,250,32769) # Spected range of the antenna, edit if needed.
eta_nu = eta(Freqs) # There should be a file named eta_nu.dat in the directory with the efficiency of the antenna.
mask = (Freqs>=50)&(Freqs<=91) # Edit if other bandrange is needed.
bwidth=1. # In MHz.

def Filter(PATH,Output_folder='.',outcome=0.):
    """
    Reduces the data to the selected band range (50-90 MHz) calculating the temperature of the measurements
    and stores them in HDF5 files. Data should be a table of the masurements in dBm's.
        
    It is spected that the files are stores in folders named with the hours of the data
    measurement.
    Example: 2013-06-13-23, 2013-06-14-00, 2013-06-14-01.
    
    Also the files inside these folders should end with .dat extension with different 
    names for every measure type, i.e short.dat, 50ohm.dat, antenna.dat, etc.
    By default this function ignores empty .dat files.
    
    Parameters:
    PATH: Path of the directory with the data files.
    Output_folfer: Path of the directory where the files should be stored.

    Optional parameters:
    outcome: Expected noise that should be extracted from the measurement, default is 0.
    
    Output:
    HDF5 tables of the temperature for the corresponding measure of the folders given.
    Also stores the gain K for the Johnson-noise calibration method.
    Also stores the ambient temperature given in the header of the antenna.dat files, 
    if the header has no ambient temperature this script will stop. 
    It is assumed that the header of the file is 19 lines long and that the last line of
    this header has the temperature in Celsius, if this is not the case modify the Tamb line.
    
    
    """
    folders = glob.glob(PATH+'/*')
    folders.sort()
    i=-1
    
    # Create target directories
    if not os.path.exists(Output_folder+'/short'):
        os.makedirs(Output_folder+'/short')
    if not os.path.exists(Output_folder+'/50ohm'):
        os.makedirs(Output_folder+'/50ohm')  
    if not os.path.exists(Output_folder+'/antenna'):
        os.makedirs(Output_folder+'/antenna')
    if not os.path.exists(Output_folder+'/Tmeas'):
        os.makedirs(Output_folder+'/Tmeas')   
    if not os.path.exists(Output_folder+'/K_jnc'): 
        os.makedirs(Output_folder+'/K_jnc')
    
    for subdirs, dirs, files in os.walk(PATH):
        dirs[:] = [d for d in dirs if not d.startswith('.')] # Inore hidden folders (ipynb checkpoints for example)
        dirs.sort()
        files.sort()
        short,antenna,_50ohm,measure,K_jnc = [],[],[],[],[]
        short_date,_50ohm_date,measure_date =[],[],[]

        # Walk through directories
        for file in files:
            path = os.path.join(subdirs,file)
            date = file.split("_")[0]
            if  os.path.getsize(path)==0: # Filtering empty data
                print 'EMPTY FILE:',path
                continue
                
            data =  np.loadtxt(path,unpack=True)
            if data.size == 0:
                print 'NO DATA IN FILE:',path
                continue
                
            elif file.endswith('short.dat'):
                T_short = Res2Temp(data,bwidth)
                short.append(T_short),short_date.append(date)
            elif file.endswith('50ohm.dat'):
                T_50ohm = Res2Temp(data,bwidth)
                _50ohm.append(T_50ohm),_50ohm_date.append(date)
            elif file.endswith('noise.dat'):
                dB_noise = data
            elif file.endswith('antenna.dat'):
                dB_antenna = data
                dB_clean = dB_antenna - dB_noise - outcome
                T_antenna = Radio_source_trans(dB_clean, Freqs, bwidth)
                T_measure = T_antenna/eta_nu - T_short # Uncalibrated measure
                Tamb = round(np.genfromtxt(path,comments='!',skip_header= 18,max_rows=1)[1]+273.15,2)
                Kjnc = Tamb/(T_50ohm-T_short) # Johnson-noise calibration coefficient
                antenna.append(T_antenna),measure.append(T_measure),K_jnc.append(Kjnc)
                measure_date.append(date)
            
        # HDF5 Table Generation       
        if i>=0 and i<len(folders) and short and antenna and _50ohm and measure and K_jnc:
            name = os.path.normpath(folders[i])
            name = name.split("/")[1]
            short = np.transpose(short)
            antenna = np.transpose(antenna)
            _50ohm = np.transpose(_50ohm)
            measure = np.transpose(measure)
            K_jnc = np.transpose(K_jnc)

            short_table = pd.DataFrame(short[mask], index = Freqs[mask], columns = short_date)
            short_table.to_hdf(Output_folder+'/short/'+name+'.hdf5','df')
            _50ohm_table = pd.DataFrame(_50ohm[mask], index = Freqs[mask], columns = _50ohm_date)
            _50ohm_table.to_hdf(Output_folder+'/50ohm/'+name+'.hdf5','df')
            antenna_table = pd.DataFrame(antenna[mask], index = Freqs[mask], columns = measure_date)
            antenna_table.to_hdf(Output_folder+'/antenna/'+name+'.hdf5','df')
            measure_table = pd.DataFrame(measure[mask], index = Freqs[mask], columns = measure_date)
            measure_table.to_hdf(Output_folder+'/Tmeas/'+name+'.hdf5','df')
            Kjnc_table = pd.DataFrame(K_jnc[mask], index = Freqs[mask], columns = measure_date)
            Kjnc_table.to_hdf(Output_folder+'/K_jnc/'+name+'.hdf5','df')
        i+=1