import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

#====================
diff_LWD=np.load('diff_LWD.npy')
diff_LWU=np.load('diff_LWU.npy')
diff_SWD=np.load('diff_SWD.npy')
diff_SWU=np.load('diff_SWU.npy')
diff_SHF=np.load('diff_SHF.npy')
diff_LHF=np.load('diff_LHF.npy')

#====================
fig, ax = plt.subplots()

ax.plot(np.arange(14), diff_LWD,'-',label='Downward Longwave',linewidth=1.0,color='tab:brown',zorder=2)
ax.plot(np.arange(14),-diff_LWU,'--',label='Upward Longwave',linewidth=1.0,color='tab:brown',zorder=3)
ax.plot(np.arange(14), diff_SWD,'-',label='Downward Shortwave',linewidth=1.0,color='tab:orange',zorder=4)
ax.plot(np.arange(14),-diff_SWU,'--',label='Upward Shortwave',linewidth=1.0,color='tab:orange',zorder=5)
ax.plot(np.arange(14), diff_SHF,'-',label='Sensible Heat Flux',linewidth=1.0,color='tab:cyan',zorder=6)
ax.plot(np.arange(14), diff_LHF,'-',label='Latent Heat Flux',linewidth=1.0,color='darkblue',zorder=7)
ax.set_xticks(np.arange(1,13,1))
ax.set_xticklabels(['M','J','J','A','S','O','N','D','J','F','M','A'])
plt.xlim(0.5,12.5)
ax.legend(fontsize=8)
ax.plot([0.5,12.5],[0,0],'--',linewidth=0.5,color='black',zorder=1)
ax.set_ylabel('Energy flux anomaly (W.m-2)',size=12)
ax.tick_params(axis='both', which='major', labelsize=10)

fig.savefig('diff_energy_seasonal_cycle_ice_shelves_ERAi_CMIP5.pdf')
