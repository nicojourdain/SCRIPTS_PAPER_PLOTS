import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import xarray as xr
import os


file_z='maps_proj_2100.npz'

fig, axs = plt.subplots(nrows=2,ncols=2,figsize=(18.0,18.0))
axs = axs.ravel()
plt.subplots_adjust(left=0.01,right=0.99,hspace=0.02,bottom=0.01,top=0.99,wspace=0.02)
npt=60 # nb of points removed in a halo (to make map appear larger)
npt2=30

# MAR masks :
grd=xr.open_dataset('RCM_ice_regrid_04000m.nc2',decode_cf=False)
msk_isf = ( grd.ICE_MAR - grd.GROUND ) * grd.af2
msk_grd = grd.GROUND * grd.af2
msk_ice = grd.AIS.where( (grd.AIS>0.) )* 0 + 1
print(msk_ice.max(),msk_ice.min())

# Topo (just for the plot):
topo=xr.open_dataset('/data/njourdain/DATA_ISMIP6/bedmap2_8km.nc')
topo2=xr.open_dataset('/data/njourdain/DATA_ISMIP6/BedMachineAntarctica_2020-07-15_v02_8km.nc')

model  = [ 'ACCESS-CM2', 'ACCESS-ESM1-5', 'CanESM5' , 'CESM2'    , 'CESM2-WACCM', 'CNRM-CM6-1', 'CNRM-ESM2-1', 'GFDL-CM4', 'GFDL-ESM4', 'GISS-E2-1-H', 'INM-CM5-0', 'IPSL-CM6A-LR', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'NorESM2-MM', 'UKESM1-0-LL' ]
member = [ 'r1i1p1f1'  , 'r1i1p1f1'     , 'r1i1p1f1', 'r11i1p1f1', 'r1i1p1f1'   , 'r1i1p1f2'  , 'r1i1p1f2'   , 'r1i1p1f1', 'r1i1p1f1' , 'r1i1p1f2'   , 'r1i1p1f1' , 'r1i1p1f1'    , 'r1i1p1f1'     , 'r1i1p1f1'  , 'r1i1p1f1'  , 'r1i1p1f2'    ]
weight = [      11     ,      24        ,       3   ,        6   ,      10      ,      10     ,      10      ,      24   ,      47    ,      41      ,      18    ,      12       ,      43        ,      39    ,       47     ,       5       ]
Nmod = np.size(model)

map_1851_1870 = np.zeros(np.shape(msk_grd))
map_2081_2100_ssp126 = np.zeros(np.shape(msk_grd))
map_2081_2100_ssp245 = np.zeros(np.shape(msk_grd))
map_2081_2100_ssp585 = np.zeros(np.shape(msk_grd))
wgt_1851_1870 = 0
wgt_2081_2100_ssp126 = 0
wgt_2081_2100_ssp245 = 0
wgt_2081_2100_ssp585 = 0

############################################################

if not os.path.exists(file_z): 

   for kmod in np.arange(Nmod):
   
      print(model[kmod])
   
      file_1850 = 'MAR-'+model[kmod]+'_asmb_1850-1979_histo_regrid_04000m_EXTENDED.nc'
      if not os.path.exists(file_1850):
        file_1850 = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_1850-1979_histo_regrid_04000m_FROM_6_MODELS_EXTENDED.nc'
   
      file_his = 'MAR-'+model[kmod]+'_asmb_1980-2014_histo_regrid_04000m.nc'
      if not os.path.exists(file_his):
        file_his = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_1980-2014_histo_regrid_04000m_FROM_6_MODELS.nc'
   
      file_ssp126 = 'MAR-'+model[kmod]+'_asmb_2015-2100_ssp126_regrid_04000m.nc'
      if not os.path.exists(file_ssp126):
        file_ssp126 = 'MAR-'+model[kmod]+'_asmb_2015-2100_ssp126_regrid_04000m_FROM_ssp585.nc'
        if not os.path.exists(file_ssp126):
          file_ssp126 = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_2015-2100_ssp126_regrid_04000m_FROM_6_MODELS.nc'
   
      file_ssp245 = 'MAR-'+model[kmod]+'_asmb_2015-2100_ssp245_regrid_04000m.nc'
      if not os.path.exists(file_ssp245):
        file_ssp245 = 'MAR-'+model[kmod]+'_asmb_2015-2100_ssp245_regrid_04000m_FROM_ssp585.nc'
        if not os.path.exists(file_ssp245):
          file_ssp245 = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_2015-2100_ssp245_regrid_04000m_FROM_6_MODELS.nc'
   
      file_ssp585 = 'MAR-'+model[kmod]+'_asmb_2015-2100_ssp585_regrid_04000m.nc'
      if not os.path.exists(file_ssp585):
        file_ssp585 = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_2015-2100_ssp585_regrid_04000m_FROM_6_MODELS.nc'
   
      d1=xr.open_dataset(file_1850,decode_cf=False)
      print(d1.asmb.isel(time=slice(1,21)).shape)
      map_1851_1870 = map_1851_1870 + d1.asmb.isel(time=slice(1,21)).mean(dim=["time"]).values * weight[kmod]
      wgt_1851_1870 = wgt_1851_1870 + weight[kmod]  
    
      if os.path.exists(file_ssp126): # all but GFDL-CM4 (no ssp126)
        d3=xr.open_dataset(file_ssp126,decode_cf=False)
        map_2081_2100_ssp126 = map_2081_2100_ssp126 + d3.asmb.isel(time=slice(66,86)).mean(dim=["time"]).values * weight[kmod]
        wgt_2081_2100_ssp126 = wgt_2081_2100_ssp126 + weight[kmod]   
   
      d4=xr.open_dataset(file_ssp245,decode_cf=False)
      print(d4.asmb.isel(time=slice(66,86)).shape)
      map_2081_2100_ssp245 = map_2081_2100_ssp245 + d4.asmb.isel(time=slice(66,86)).mean(dim=["time"]).values * weight[kmod]
      wgt_2081_2100_ssp245 = wgt_2081_2100_ssp245 + weight[kmod]  
    
      d5=xr.open_dataset(file_ssp585,decode_cf=False)
      map_2081_2100_ssp585 = map_2081_2100_ssp585 + d5.asmb.isel(time=slice(66,86)).mean(dim=["time"]).values * weight[kmod]
      wgt_2081_2100_ssp585 = wgt_2081_2100_ssp585 + weight[kmod] 
   
   map_1851_1870 = map_1851_1870 / wgt_1851_1870  
   map_2081_2100_ssp126 = map_2081_2100_ssp126 / wgt_2081_2100_ssp126
   map_2081_2100_ssp245 = map_2081_2100_ssp245 / wgt_2081_2100_ssp245
   map_2081_2100_ssp585 = map_2081_2100_ssp585 / wgt_2081_2100_ssp585
   print(wgt_1851_1870, wgt_2081_2100_ssp126, wgt_2081_2100_ssp245, wgt_2081_2100_ssp585)

   np.savez(file_z,map_1851_1870 = map_1851_1870, map_2081_2100_ssp126 = map_2081_2100_ssp126, map_2081_2100_ssp245 = map_2081_2100_ssp245, map_2081_2100_ssp585 = map_2081_2100_ssp585)

else:

   zz=np.load(file_z)
   map_1851_1870 = zz['map_1851_1870']
   map_2081_2100_ssp126 = zz['map_2081_2100_ssp126']
   map_2081_2100_ssp245 = zz['map_2081_2100_ssp245']
   map_2081_2100_ssp585 = zz['map_2081_2100_ssp585']

map_1851_1870 = map_1851_1870 * msk_ice.values * 365.25 * 86400      # kg m-2 yr-1
map_2081_2100_ssp126 = map_2081_2100_ssp126 * msk_ice.values * 365.25 * 86400
map_2081_2100_ssp245 = map_2081_2100_ssp245 * msk_ice.values * 365.25 * 86400
map_2081_2100_ssp585 = map_2081_2100_ssp585 * msk_ice.values * 365.25 * 86400
 
##########################################################################
# PLOT :

#----------
# Colorbar:
cbar_range = np.arange(-275.,300.,25.)

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

cax=fig.add_axes([0.35, 0.53, 0.45, 0.02]) # color bar
im0=axs[0].contourf(grd.x[npt:-npt],grd.y[npt:-npt],map_1851_1870[npt:-npt,npt:-npt],cbar_range,cmap=cmap_new,extend='both')
cbar0=fig.colorbar(im0, cax=cax, orientation="horizontal")
cbar0.ax.tick_params(labelsize=16)
cbar0.set_label(r'kg m$^{-2}$ yr$^{-1}$',labelpad=-70,fontsize=16)
axs[0].contour(topo.x[npt2:-npt2],topo.y[npt2:-npt2],topo.surface.values[npt2:-npt2,npt2:-npt2],[1000,2000,3000,4000],colors='darkgrey',linewidths=1)
axs[0].contour(topo2.x[npt2:-npt2],topo2.y[npt2:-npt2],topo2.mask.values[npt2:-npt2,npt2:-npt2],[0.5],colors='black',linewidths=1.5)
xc=np.mean(axs[0].get_xlim())
yu=axs[0].get_ylim()[1]
axs[0].text(xc,0.95*yu,'(a) 1851-1870 SMB anomaly w.r.t. 1995-2014',fontsize=18,fontweight='bold',ha='center')
axs[0].set_axis_off()

im1=axs[1].contourf(grd.x[npt:-npt],grd.y[npt:-npt],map_2081_2100_ssp126[npt:-npt,npt:-npt],cbar_range,cmap=cmap_new,extend='both')
axs[1].contour(topo.x[npt2:-npt2],topo.y[npt2:-npt2],topo.surface.values[npt2:-npt2,npt2:-npt2],[1000,2000,3000,4000],colors='darkgrey',linewidths=1)
axs[1].contour(topo2.x[npt2:-npt2],topo2.y[npt2:-npt2],topo2.mask.values[npt2:-npt2,npt2:-npt2],[0.5],colors='black',linewidths=1.5)
axs[1].text(xc,0.95*yu,'(b) SSP1-2.6 2081-2100 SMB anomaly w.r.t. 1995-2014',fontsize=18,fontweight='bold',ha='center')
axs[1].set_axis_off()

im2=axs[2].contourf(grd.x[npt:-npt],grd.y[npt:-npt],map_2081_2100_ssp245[npt:-npt,npt:-npt],cbar_range,cmap=cmap_new,extend='both')
axs[2].contour(topo.x[npt2:-npt2],topo.y[npt2:-npt2],topo.surface.values[npt2:-npt2,npt2:-npt2],[1000,2000,3000,4000],colors='darkgrey',linewidths=1)
axs[2].contour(topo2.x[npt2:-npt2],topo2.y[npt2:-npt2],topo2.mask.values[npt2:-npt2,npt2:-npt2],[0.5],colors='black',linewidths=1.5)
axs[2].text(xc,0.95*yu,'(c) SSP2-4.5 2081-2100 SMB anomaly w.r.t. 1995-2014',fontsize=18,fontweight='bold',ha='center')
axs[2].set_axis_off()

im3=axs[3].contourf(grd.x[npt:-npt],grd.y[npt:-npt],map_2081_2100_ssp585[npt:-npt,npt:-npt],cbar_range,cmap=cmap_new,extend='both')
axs[3].contour(topo.x[npt2:-npt2],topo.y[npt2:-npt2],topo.surface.values[npt2:-npt2,npt2:-npt2],[1000,2000,3000,4000],colors='darkgrey',linewidths=1)
axs[3].contour(topo2.x[npt2:-npt2],topo2.y[npt2:-npt2],topo2.mask.values[npt2:-npt2,npt2:-npt2],[0.5],colors='black',linewidths=1.5)
axs[3].text(xc,0.95*yu,'(d) SSP5-8.5 2081-2100 SMB anomaly w.r.t. 1995-2014',fontsize=18,fontweight='bold',ha='center')
axs[3].set_axis_off()

##########################################################################

fig.savefig("map_projection_16_models.pdf")
