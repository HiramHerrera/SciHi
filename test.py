import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from power_to_temperature import *

f_eta, p_eta = np.loadtxt('eta_nu.dat', unpack=True)
f_eta *= 1e-6
p_eta *= 0.01

def bandwidth(min_freq, max_freq):
    """
    Funcion para encontrar el binning en las frecuencias para los datos de
    sci-hi
    """

    n = 32769.0
    top_freq = 250.0 # MHz
    val = n / top_freq

    low = round(min_freq * val)
    top = round(max_freq * val)


    low_val = low * val**-1
    top_val = top * val**-1
    bin = top - low

    n = 32769.0
    top_freq = 250.0 # MHz
    val = n / top_freq
    Bwidth = (top_val - low_val) / bin #MHz
    Bwidth *= 1e6 #Hz

    fl = int(low_val * val)
    ft = int(top_val * val)

    values = np.array([fl, ft])

    frequencies = np.linspace(low_val, top_val, bin)
    return frequencies, values, Bwidth


def eta(frequencies):
    """
    Interpolation of the efficiency of the antenna to match our bandwith
    """
    z = interpolate.splrep(f_eta, p_eta, s=0)
    efficiency = interpolate.splev(frequencies, z, der = 0)
    return efficiency

def calibration(antenna_data, noise_data, short_data, res_data, temp, freqs, Bwidth):
    """
    calibration function. Voytek et al. (2014)
    """

    antenna_power = antenna_data  - noise_data
    T_antenna = Radio_source_trans(antenna_power, freqs, Bwidth) #K
    T_short = Res2Temp(short_data, Bwidth) #K
    T_res = Res2Temp(res_data, Bwidth)  #K

#mean_sky_temp = np.mean(T_antenna, axis=0)
    K = temp / (T_res - T_short)

    power = T_antenna/ eta(freqs)
    Temperature = K * (power - T_short)

    return Temperature
