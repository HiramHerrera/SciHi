import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

PATH='/Users/hiram/Codes/21cmFAST-master/Output_files/Ts_outs/'
file = PATH+'global_evolution_zetaIon30.00_Nsteps40_zprimestepfactor1.020_zetaX2.0e+56_alphaX1.2_TvirminX3.0e+04_Pop2_32_300Mpc'


data = open(file,'r')
lines = data.readlines()
data.close()

dat=[]
for line in lines:
	dat.append(line.split())
dat = np.array(dat,dtype=np.float) 
dat = np.asfarray(dat,float)
z = dat[:,0]
Xh = dat[:,9]
Ts = dat[:,4]
Tgamma = dat[:,5]/(1.+z)

Tb = 27.*(Xh)*np.sqrt((1.+z)/10.)*(1.- Tgamma/Ts)

plt.plot(z,Ts,color='k',label=r'$T_s$')
plt.plot(z,Tgamma,color='b',label=r'$T_\gamma$')
plt.legend()
plt.savefig('zvsTsTgam.png')
plt.close()

plt.plot(z,Tb)
plt.savefig('zvsdTb.png')

