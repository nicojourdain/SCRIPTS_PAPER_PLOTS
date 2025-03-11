import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

fig, axs = plt.subplots(nrows=1,ncols=2,figsize=(18.0,9.0))
axs = axs.ravel()
plt.subplots_adjust(left=0.01,right=0.99,hspace=0.02,bottom=0.01,top=0.99,wspace=0.02)
npt=60 # nb of points removed in a halo (to make map appear larger)
npt2=30

# input files
grd=xr.open_dataset('RCM_ice_regrid_04000m.nc2',decode_cf=False)
msk_ice = grd.AIS.where( (grd.AIS>0.) )* 0 + 1
#-
sm585=xr.open_dataset('MAR-IPSL-CM6A-LR_asmb_2015-2200_ssp585_regrid_04000m.nc',decode_cf=False)
smEXT060=xr.open_dataset('MAR-IPSL-CM6A-LR_asmb_2101-2300_ssp585_regrid_04000m_EXTENDED_MOA060.nc',decode_cf=False)
#- Topo (just for the plot):
topo=xr.open_dataset('/data/njourdain/DATA_ISMIP6/bedmap2_8km.nc')
topo2=xr.open_dataset('/data/njourdain/DATA_ISMIP6/BedMachineAntarctica_2020-07-15_v02_8km.nc')

smb_MAR = sm585.asmb.isel(time=slice(85,104)).mean(dim="time") * msk_ice * 365.25 * 86400
smb_EXT = smEXT060.asmb.isel(time=slice(0,19)).mean(dim="time") * msk_ice * 365.25 * 86400

# moving the zero of colorbar
# NB: modify the Ncool to Nwarm ratio (total=256) to place zero as desired.
cbar_range=np.arange(-500,550,50)
Ncool=int(256*(-np.amin(cbar_range)/(np.amax(cbar_range)-np.amin(cbar_range))))
Nwarm=256-Ncool
col = cm.get_cmap('PuOr', 256)
tmp1 = col(np.linspace(0.47, 1.00, Ncool)) # decrease first number to have more white in the middle light-blue colors
tmp2 = col(np.linspace(0.00, 0.51, Nwarm)) # increase second number to have more white in the middle light-yellow colors
newcolors = np.append(tmp1[::-1,:],tmp2[::-1,:],axis=0)
cmap_new = ListedColormap(newcolors)

#----- a -----
cax=fig.add_axes([0.2, 0.05, 0.6, 0.03]) # color bar
im0=axs[0].contourf(grd.x[npt:-npt],grd.y[npt:-npt],smb_MAR[npt:-npt,npt:-npt],cbar_range,cmap=cmap_new,extend='both')
cbar0=fig.colorbar(im0, cax=cax, orientation="horizontal")
cbar0.ax.tick_params(labelsize=13)
cbar0.set_label(r'kg m$^{-2}$ yr$^{-1}$',labelpad=-60,fontsize=14)
axs[0].contour(topo.x[npt2:-npt2],topo.y[npt2:-npt2],topo.surface.values[npt2:-npt2,npt2:-npt2],[1000,2000,3000,4000],colors='darkgrey',linewidths=1)
axs[0].contour(topo2.x[npt2:-npt2],topo2.y[npt2:-npt2],topo2.mask.values[npt2:-npt2,npt2:-npt2],[0.5],colors='black',linewidths=1.5)
xc=np.mean(axs[0].get_xlim())
yu=axs[0].get_ylim()[1]
axs[0].text(xc,0.92*yu,'(a) SMB anomaly in the original MAR simulation',fontsize=16,fontweight='bold',ha='center')
axs[0].set_axis_off()

#----- b -----
im1=axs[1].contourf(grd.x[npt:-npt],grd.y[npt:-npt],smb_EXT[npt:-npt,npt:-npt],cbar_range,cmap=cmap_new,extend='both')
axs[1].contour(topo.x[npt2:-npt2],topo.y[npt2:-npt2],topo.surface.values[npt2:-npt2,npt2:-npt2],[1000,2000,3000,4000],colors='darkgrey',linewidths=1)
axs[1].contour(topo2.x[npt2:-npt2],topo2.y[npt2:-npt2],topo2.mask.values[npt2:-npt2,npt2:-npt2],[0.5],colors='black',linewidths=1.5)
xc=np.mean(axs[1].get_xlim())
yl=axs[1].get_ylim()[1]
axs[1].text(xc,0.92*yl,'(b) emulated SMB anomaly',fontsize=16,fontweight='bold',ha='center')
axs[1].set_axis_off()

#-----

figname='map_eval_1980_to_2200.pdf'

fig.savefig(figname)

