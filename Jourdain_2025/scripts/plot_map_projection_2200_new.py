import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import xarray as xr
import os


file_z='maps_proj_2200.npz'

fig, axs = plt.subplots(nrows=2,ncols=2,figsize=(18.0,18.0))
axs = axs.ravel()
plt.subplots_adjust(left=0.01,right=0.99,hspace=0.02,bottom=0.01,top=0.99,wspace=0.02)
npt=60 # nb of points removed in a halo (to make map appear larger)
npt2=30

# masks:
grd=xr.open_dataset('RCM_ice_regrid_04000m.nc2',decode_cf=False)
msk_isf = ( grd.ICE_MAR - grd.GROUND ) * grd.af2
msk_grd = grd.GROUND * grd.af2
msk_ice = grd.AIS.where( (grd.AIS>0.) )* 0 + 1

# Topo (just for the plot):
topo=xr.open_dataset('/data/njourdain/DATA_ISMIP6/bedmap2_8km.nc')
topo2=xr.open_dataset('/data/njourdain/DATA_ISMIP6/BedMachineAntarctica_2020-07-15_v02_8km.nc')

model  = [ 'ACCESS-CM2', 'ACCESS-ESM1-5', 'CanESM5' , 'CESM2-WACCM', 'GISS-E2-1-H', 'IPSL-CM6A-LR',  'MRI-ESM2-0', 'UKESM1-0-LL' ]
member = [ 'r1i1p1f1'  , 'r1i1p1f1'     , 'r1i1p1f1', 'r1i1p1f1'   , 'r1i1p1f2'   , 'r1i1p1f1'    ,  'r1i1p1f1'  , 'r4i1p1f2'    ]
Nmod = np.size(model)

map_2181_2200_A = np.zeros(np.shape(msk_grd))
map_2181_2200_B = np.zeros(np.shape(msk_grd))
nA = 0
nB = 0

############################################################

if not os.path.exists(file_z): 


   for kmod in np.arange(Nmod):
      
     print(model[kmod])
        
     file_his = 'MAR-'+model[kmod]+'_asmb_1980-2014_histo_regrid_04000m.nc'
     if not os.path.exists(file_his):
       file_his = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_1980-2014_histo_regrid_04000m_FROM_6_MODELS.nc'
     if ( model[kmod] == 'UKESM1-0-LL' ):
       file_his = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_1980-2014_histo_regrid_04000m_FROM_UKESM1-0-LL-r1i1p1f2-histo.nc'
      
     file_ssp126 = 'MAR-'+model[kmod]+'_asmb_2015-2200_ssp126_regrid_04000m.nc'
     if not os.path.exists(file_ssp126):
       file_ssp126 = 'MAR-'+model[kmod]+'_asmb_2015-2200_ssp126_regrid_04000m_FROM_ssp585.nc'
       if not os.path.exists(file_ssp126):
         file_ssp126 = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_2015-2200_ssp126_regrid_04000m_MERGED.nc'
        
     file_ssp585 = 'MAR-'+model[kmod]+'_asmb_2015-2200_ssp585_regrid_04000m.nc'
     if not os.path.exists(file_ssp585):
       file_ssp585 = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_2015-2200_ssp585_regrid_04000m_MERGED.nc'
        
     dd=xr.open_dataset(file_ssp585,decode_cf=False)
     if ( ( model[kmod] == 'CanESM5' ) | ( model[kmod] == 'CESM2-WACCM' ) | ( model[kmod] == 'IPSL-CM6A-LR' ) | ( model[kmod] == 'UKESM1-0-LL' ) ):
       map_2181_2200_A = map_2181_2200_A + dd.asmb.isel(time=slice(165,186)).mean(dim=["time"]).values
       nA = nA + 1
     else:
       map_2181_2200_B = map_2181_2200_B + dd.asmb.isel(time=slice(165,186)).mean(dim=["time"]).values
       nB = nB + 1
  
   print(nA, nB)
   map_2181_2200_A = map_2181_2200_A / nA
   map_2181_2200_B = map_2181_2200_B / nB

   np.savez(file_z,map_2181_2200_A = map_2181_2200_A, map_2181_2200_B = map_2181_2200_B)

else:

   zz=np.load(file_z)
   map_2181_2200_A = zz['map_2181_2200_A']
   map_2181_2200_B = zz['map_2181_2200_B']

map_2181_2200_A = map_2181_2200_A * msk_ice.values * 365.25 * 86400      # mm w. eq. / year
map_2181_2200_B = map_2181_2200_B * msk_ice.values * 365.25 * 86400      # mm w. eq. / year

##########################################################################
# Temperature anomalies

b1=xr.open_dataset('/data/njourdain/DATA_PROTECT/TAS/tas_Ayr_ACCESS-CM2_ssp585_r1i1p1f1_2101_2300.nc',decode_cf=False)
b2=xr.open_dataset('/data/njourdain/DATA_PROTECT/TAS/tas_Ayr_ACCESS-ESM1-5_ssp585_r1i1p1f1_2101_2300.nc',decode_cf=False)
a1=xr.open_dataset('/data/njourdain/DATA_PROTECT/TAS/tas_Ayr_CanESM5_ssp585_r1i1p1f1_2101_2300.nc',decode_cf=False)
a2=xr.open_dataset('/data/njourdain/DATA_PROTECT/TAS/tas_Ayr_CESM2-WACCM_ssp585_r1i1p1f1_2101_2300.nc',decode_cf=False)
b3=xr.open_dataset('/data/njourdain/DATA_PROTECT/TAS/tas_Ayr_GISS-E2-1-H_ssp585_r1i1p1f2_2101_2300.nc',decode_cf=False)
a3=xr.open_dataset('/data/njourdain/DATA_PROTECT/TAS/tas_Ayr_IPSL-CM6A-LR_ssp585_r1i1p1f1_2101_2300.nc',decode_cf=False)
b4=xr.open_dataset('/data/njourdain/DATA_PROTECT/TAS/tas_Ayr_MRI-ESM2-0_ssp585_r1i1p1f1_2101_2300.nc',decode_cf=False)
a4=xr.open_dataset('/data/njourdain/DATA_PROTECT/TAS/tas_Ayr_UKESM1-0-LL_ssp585_r4i1p1f2_2101_2300.nc',decode_cf=False)

climb1=xr.open_dataset('/data/njourdain/DATA_PROTECT/TAS/tas_Ayr_ACCESS-CM2_historical_r1i1p1f1_1850_2014.nc',decode_cf=False)
climb2=xr.open_dataset('/data/njourdain/DATA_PROTECT/TAS/tas_Ayr_ACCESS-ESM1-5_historical_r1i1p1f1_1850_2014.nc',decode_cf=False)
clima1=xr.open_dataset('/data/njourdain/DATA_PROTECT/TAS/tas_Ayr_CanESM5_historical_r1i1p1f1_1850_2014.nc',decode_cf=False)
clima2=xr.open_dataset('/data/njourdain/DATA_PROTECT/TAS/tas_Ayr_CESM2-WACCM_historical_r1i1p1f1_1850_2014.nc',decode_cf=False)
climb3=xr.open_dataset('/data/njourdain/DATA_PROTECT/TAS/tas_Ayr_GISS-E2-1-H_historical_r1i1p1f2_1850_2014.nc',decode_cf=False)
clima3=xr.open_dataset('/data/njourdain/DATA_PROTECT/TAS/tas_Ayr_IPSL-CM6A-LR_historical_r1i1p1f1_1850_2014.nc',decode_cf=False)
climb4=xr.open_dataset('/data/njourdain/DATA_PROTECT/TAS/tas_Ayr_MRI-ESM2-0_historical_r1i1p1f1_1850_2014.nc',decode_cf=False)
clima4=xr.open_dataset('/data/njourdain/DATA_PROTECT/TAS/tas_Ayr_UKESM1-0-LL_historical_r4i1p1f2_1850_2014.nc',decode_cf=False)

dTa = (   a1.tas.isel(time=slice(81,100)).mean(dim=["time"]) - clima1.tas.isel(time=slice(146,165)).mean(dim=["time"])
        + a2.tas.isel(time=slice(81,100)).mean(dim=["time"]) - clima2.tas.isel(time=slice(146,165)).mean(dim=["time"])
        + a3.tas.isel(time=slice(81,100)).mean(dim=["time"]) - clima3.tas.isel(time=slice(146,165)).mean(dim=["time"])
        + a4.tas.isel(time=slice(81,100)).mean(dim=["time"]) - clima4.tas.isel(time=slice(146,165)).mean(dim=["time"]) ).values * msk_ice.values * 0.250000 

dTb = (   b1.tas.isel(time=slice(81,100)).mean(dim=["time"]) - climb1.tas.isel(time=slice(146,165)).mean(dim=["time"])
        + b2.tas.isel(time=slice(81,100)).mean(dim=["time"]) - climb2.tas.isel(time=slice(146,165)).mean(dim=["time"])
        + b3.tas.isel(time=slice(81,100)).mean(dim=["time"]) - climb3.tas.isel(time=slice(146,165)).mean(dim=["time"])
        + b4.tas.isel(time=slice(81,100)).mean(dim=["time"]) - climb4.tas.isel(time=slice(146,165)).mean(dim=["time"]) ).values * msk_ice.values * 0.250000 

print(np.shape(dTa))
print(np.shape(dTb))

##########################################################################
# PLOT :

#----------
# Colorbar:
cbar_range = np.arange(-3000.,1750.,250.)

#----------
# Defining colormap:

# moving the zero of colorbar
# NB: modify the Ncool to Nwarm ratio (total=256) to place zero as desired.
Ncool=int(256*(-np.amin(cbar_range)/(np.amax(cbar_range)-np.amin(cbar_range))))
Nwarm=256-Ncool
col = cm.get_cmap('PuOr', 256)
tmp1 = col(np.linspace(0.47, 1.00, Ncool)) # decrease first number to have more white in the middle light-blue colors
tmp2 = col(np.linspace(0.00, 0.51, Nwarm)) # increase second number to have more white in the middle light-yellow colors
newcolors = np.append(tmp1[::-1,:],tmp2[::-1,:],axis=0)
cmap_new = ListedColormap(newcolors)

cay=fig.add_axes([0.05, 0.02, 0.4, 0.015]) # color bar
im8=axs[0].contourf(grd.x.values[npt:-npt],grd.y.values[npt:-npt],dTa[npt:-npt,npt:-npt],np.arange(7,19.0,1.0),cmap='plasma',extend='both')
cbar8=fig.colorbar(im8, cax=cay, orientation="horizontal")
cbar8.ax.tick_params(labelsize=13)
cbar8.set_label('°C',labelpad=-60,fontsize=14)
axs[0].contour(topo.x[npt2:-npt2],topo.y[npt2:-npt2],topo.surface.values[npt2:-npt2,npt2:-npt2],[1000,2000,3000,4000],colors='darkgrey',linewidths=1)
axs[0].contour(topo2.x[npt2:-npt2],topo2.y[npt2:-npt2],topo2.mask.values[npt2:-npt2,npt2:-npt2],[0.5],colors='black',linewidths=1.5)
xc=np.mean(axs[0].get_xlim())
yu=axs[0].get_ylim()[1]
axs[0].text(xc,0.99*yu,'(a) SSP5-8.5 2181-2200 air temperature anomaly w.r.t. 1995-2014',fontsize=16,fontweight='bold',ha='center')
axs[0].text(xc,0.92*yu,'[CanESM5, CESM2-WACCM, IPSL-CM6A-LR, UKESM1-0-LL]',fontsize=16,fontweight='bold',ha='center')
axs[0].set_axis_off()

im9=axs[2].contourf(grd.x.values[npt:-npt],grd.y.values[npt:-npt],dTb[npt:-npt,npt:-npt],np.arange(7,19.0,1.0),cmap='plasma',extend='both')
axs[2].contour(topo.x[npt2:-npt2],topo.y[npt2:-npt2],topo.surface.values[npt2:-npt2,npt2:-npt2],[1000,2000,3000,4000],colors='darkgrey',linewidths=1)
axs[2].contour(topo2.x[npt2:-npt2],topo2.y[npt2:-npt2],topo2.mask.values[npt2:-npt2,npt2:-npt2],[0.5],colors='black',linewidths=1.5)
xc=np.mean(axs[2].get_xlim())
yl=axs[2].get_ylim()[1]
axs[2].text(xc,0.99*yl,'(c) SSP5-8.5 2181-2200 air temperature anomaly w.r.t. 1995-2014',fontsize=16,fontweight='bold',ha='center')
axs[2].text(xc,0.92*yl,'[ACCESS-CM2, ACCESS-ESM1-5, GISS-E2-1-H, MRI-ESM2-0]',fontsize=16,fontweight='bold',ha='center')
axs[2].set_axis_off()

cax=fig.add_axes([0.55, 0.02, 0.4, 0.015]) # color bar
im0=axs[1].contourf(grd.x[npt:-npt],grd.y[npt:-npt],map_2181_2200_A[npt:-npt,npt:-npt],cbar_range,cmap=cmap_new,extend='both')
cbar0=fig.colorbar(im0, cax=cax, orientation="horizontal")
cbar0.ax.tick_params(labelsize=13)
cbar0.set_label(r'kg m$^{-2}$ yr$^{-1}$',labelpad=-60,fontsize=14)
axs[1].contour(topo.x[npt2:-npt2],topo.y[npt2:-npt2],topo.surface.values[npt2:-npt2,npt2:-npt2],[1000,2000,3000,4000],colors='darkgrey',linewidths=1)
axs[1].contour(topo2.x[npt2:-npt2],topo2.y[npt2:-npt2],topo2.mask.values[npt2:-npt2,npt2:-npt2],[0.5],colors='black',linewidths=1.5)
xc=np.mean(axs[1].get_xlim())
yu=axs[1].get_ylim()[1]
axs[1].text(xc,0.99*yu,'(b) SSP5-8.5 2181-2200 SMB anomaly w.r.t. 1995-2014',fontsize=16,fontweight='bold',ha='center')
axs[1].text(xc,0.92*yu,'[CanESM5, CESM2-WACCM, IPSL-CM6A-LR, UKESM1-0-LL]',fontsize=16,fontweight='bold',ha='center')
axs[1].set_axis_off()

im1=axs[3].contourf(grd.x[npt:-npt],grd.y[npt:-npt],map_2181_2200_B[npt:-npt,npt:-npt],cbar_range,cmap=cmap_new,extend='both')
axs[3].contour(topo.x[npt2:-npt2],topo.y[npt2:-npt2],topo.surface.values[npt2:-npt2,npt2:-npt2],[1000,2000,3000,4000],colors='darkgrey',linewidths=1)
axs[3].contour(topo2.x[npt2:-npt2],topo2.y[npt2:-npt2],topo2.mask.values[npt2:-npt2,npt2:-npt2],[0.5],colors='black',linewidths=1.5)
xc=np.mean(axs[3].get_xlim())
yl=axs[3].get_ylim()[1]
axs[3].text(xc,0.99*yl,'(d) SSP5-8.5 2181-2200 SMB anomaly w.r.t. 1995-2014',fontsize=16,fontweight='bold',ha='center')
axs[3].text(xc,0.92*yl,'[ACCESS-CM2, ACCESS-ESM1-5, GISS-E2-1-H, MRI-ESM2-0]',fontsize=16,fontweight='bold',ha='center')
axs[3].set_axis_off()

##########################################################################

fig.savefig("map_projection_8_models_2200.pdf")
