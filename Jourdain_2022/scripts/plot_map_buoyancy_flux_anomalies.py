import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import stereo as st
import gsw

file_SI_p_A='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MAR/climato_monthly_AMUXL12-GNJ002_BM02MAR_icemod_1989_2009.nc'
file_SI_p_B='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM03MAR/climato_monthly_AMUXL12-GNJ002_BM03MAR_icemod_1989_2009.nc'
file_SI_p_C='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM04MAR/climato_monthly_AMUXL12-GNJ002_BM04MAR_icemod_1989_2009.nc'
file_SI_f_A='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MARrcp85/climato_monthly_AMUXL12-GNJ002_BM02MARrcp85_icemod_1989_2009.nc'
file_SI_f_B='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM03MARrcBDY/climato_monthly_AMUXL12-GNJ002_BM03MARrcBDY_icemod_1989_2009.nc'
file_SI_f_C='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM04MARrcp85/climato_monthly_AMUXL12-GNJ002_BM04MARrcp85_icemod_1989_2009.nc'

ncpA = xr.open_dataset(file_SI_p_A,decode_cf=False)
ncpB = xr.open_dataset(file_SI_p_B,decode_cf=False)
ncpC = xr.open_dataset(file_SI_p_C,decode_cf=False)
ncfA = xr.open_dataset(file_SI_f_A,decode_cf=False)
ncfB = xr.open_dataset(file_SI_f_B,decode_cf=False)
ncfC = xr.open_dataset(file_SI_f_C,decode_cf=False)

file_TS_p_A='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MAR/climato_monthly_AMUXL12-GNJ002_BM02MAR_grid_T_1989_2009.nc'
file_TS_p_B='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM03MAR/climato_monthly_AMUXL12-GNJ002_BM03MAR_grid_T_1989_2009.nc'
file_TS_p_C='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM04MAR/climato_monthly_AMUXL12-GNJ002_BM04MAR_grid_T_1989_2009.nc'
file_TS_f_A='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MARrcp85/climato_monthly_AMUXL12-GNJ002_BM02MARrcp85_grid_T_1989_2009.nc'
file_TS_f_B='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM03MARrcBDY/climato_monthly_AMUXL12-GNJ002_BM03MARrcBDY_grid_T_1989_2009.nc'
file_TS_f_C='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM04MARrcp85/climato_monthly_AMUXL12-GNJ002_BM04MARrcp85_grid_T_1989_2009.nc'

tspA = xr.open_dataset(file_TS_p_A,decode_cf=False)
tspB = xr.open_dataset(file_TS_p_B,decode_cf=False)
tspC = xr.open_dataset(file_TS_p_C,decode_cf=False)
tsfA = xr.open_dataset(file_TS_f_A,decode_cf=False)
tsfB = xr.open_dataset(file_TS_f_B,decode_cf=False)
tsfC = xr.open_dataset(file_TS_f_C,decode_cf=False)

file_SBC_p_A='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MAR/climato_monthly_AMUXL12-GNJ002_BM02MAR_SBC_1989_2009.nc'
file_SBC_p_B='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM03MAR/climato_monthly_AMUXL12-GNJ002_BM03MAR_SBC_1989_2009.nc'
file_SBC_p_C='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM04MAR/climato_monthly_AMUXL12-GNJ002_BM04MAR_SBC_1989_2009.nc'
file_SBC_f_A='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MARrcp85/climato_monthly_AMUXL12-GNJ002_BM02MARrcp85_SBC_1989_2009.nc'
file_SBC_f_B='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM03MARrcBDY/climato_monthly_AMUXL12-GNJ002_BM03MARrcBDY_SBC_1989_2009.nc'
file_SBC_f_C='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM04MARrcp85/climato_monthly_AMUXL12-GNJ002_BM04MARrcp85_SBC_1989_2009.nc'

sbpA = xr.open_dataset(file_SBC_p_A,decode_cf=False)
sbpB = xr.open_dataset(file_SBC_p_B,decode_cf=False)
sbpC = xr.open_dataset(file_SBC_p_C,decode_cf=False)
sbfA = xr.open_dataset(file_SBC_f_A,decode_cf=False)
sbfB = xr.open_dataset(file_SBC_f_B,decode_cf=False)
sbfC = xr.open_dataset(file_SBC_f_C,decode_cf=False)

print(ncpA.vfxice.shape)

# Average buoyancy flux over A, B, C simulations ( kg . m^-1 . s^-3 ) :
g = 9.81
cp = 3992.

hfxout_p = ( ncpA.hfxout.mean(axis=0) + ncpB.hfxout.mean(axis=0) + ncpC.hfxout.mean(axis=0) ) /3.
hfxout_f = ( ncfA.hfxout.mean(axis=0) + ncfB.hfxout.mean(axis=0) + ncfC.hfxout.mean(axis=0) ) /3.

empmr_p = ( sbpA.empmr.mean(axis=0) + sbpB.empmr.mean(axis=0) + sbpC.empmr.mean(axis=0) ) /3.
empmr_f = ( sbfA.empmr.mean(axis=0) + sbfB.empmr.mean(axis=0) + sbfC.empmr.mean(axis=0) ) /3.

SST_p = ( tspA.toce.isel(deptht=0).mean(axis=0) + tspB.toce.isel(deptht=0).mean(axis=0) + tspC.toce.isel(deptht=0).mean(axis=0) ) / 3.
SST_f = ( tsfA.toce.isel(deptht=0).mean(axis=0) + tsfB.toce.isel(deptht=0).mean(axis=0) + tsfC.toce.isel(deptht=0).mean(axis=0) ) / 3.

SSS_p = ( tspA.soce.isel(deptht=0).mean(axis=0) + tspB.soce.isel(deptht=0).mean(axis=0) + tspC.soce.isel(deptht=0).mean(axis=0) ) / 3.
SSS_f = ( tsfA.soce.isel(deptht=0).mean(axis=0) + tsfB.soce.isel(deptht=0).mean(axis=0) + tsfC.soce.isel(deptht=0).mean(axis=0) ) / 3.

B_p = g/cp * gsw.density.alpha(SSS_p,SST_p,0.) * hfxout_p - g * gsw.density.beta(SSS_p,SST_p,0.) * SSS_p * empmr_p
B_f = g/cp * gsw.density.alpha(SSS_f,SST_f,0.) * hfxout_f - g * gsw.density.beta(SSS_f,SST_f,0.) * SSS_f * empmr_f

print('mean alpha:',gsw.density.alpha(SSS_p,SST_p,0.).mean().values)
print('mean beta:',gsw.density.beta(SSS_p,SST_p,0.).mean().values)

#file_msh = '/Users/njourdain/OUTPUT_AMU/mesh_mask_AMUXL12_BedMachineAntarctica-2020-07-15_v02.nc'
file_msh = '/Users/njourdain/OUTPUT_AMU/mesh_mask_AMUXL12_BedMachineAntarctica-2020-07-15_v02_ICB380.nc'
file_bat = '/Users/njourdain/OUTPUT_AMU/bathy_meter_AMUXL12_BedMachineAntarctica-2020-07-15_v02.nc'

ncmsh = xr.open_dataset(file_msh,decode_cf=False)
ncbat = xr.open_dataset(file_bat,decode_cf=False)

lon2d = ncmsh.glamt.isel(t=0)
lat2d = ncmsh.gphit.isel(t=0)

bathy = ncbat.Bathymetry

msk = ncmsh.tmask.isel(t=0,z=0)*1.0 + ncmsh.tmask.isel(t=0).max('z')*1.0 # 2 open ocean, 1 ice shelf, 0 grounded
msk = msk.where( msk<1.5, np.nan )

#-------------------------------------------------------------------------
# Average over the continental shelf:

msk_shelf = ncmsh.tmask.isel(t=0,z=0).where((ncmsh.tmask.isel(t=0,z=0)==1)&(bathy<1500)&(lon2d>=-135)&(lon2d<=-100.0)&(lat2d<-70.2),0.e0)
tmp_B_p = B_p * msk_shelf *  ncmsh.e1t.isel(t=0) * ncmsh.e2t.isel(t=0)
tmp_B_f = B_f * msk_shelf *  ncmsh.e1t.isel(t=0) * ncmsh.e2t.isel(t=0)
tmp_area = msk_shelf *  ncmsh.e1t.isel(t=0) * ncmsh.e2t.isel(t=0)

mean_B_p = tmp_B_p.sum() / tmp_area.sum()
mean_B_f = tmp_B_f.sum() / tmp_area.sum()

#-------------------------------------------------------------------------

proj=ccrs.SouthPolarStereo(central_longitude=-115.0)
trans=ccrs.PlateCarree()

fig, axs = plt.subplots(nrows=1,ncols=2,
                        subplot_kw={'projection': proj},
                        figsize=(21.0,7.0))
axs = axs.ravel()

#--------------------

# (lonmin, lonmax, latmin, latmax) :
axs[0].set_extent( (-145, -85, -76, -68), trans )

# Add and customize meridians/parallels
gl0 = axs[0].gridlines(draw_labels=True, linewidth=0.75, color='gray', alpha=0.7, linestyle='--',
                       dms=True, x_inline=False, y_inline=False )
gl0.top_labels = False
gl0.right_labels = False
gl0.xlocator = mticker.FixedLocator([-140, -130, -120, -110, -100, -90])
gl0.ylocator = mticker.FixedLocator([-78, -76, -74, -72, -70, -68])
gl0.xformatter = LongitudeFormatter()
gl0.yformatter = LatitudeFormatter()
gl0.xlabel_style = {'size': 12, 'color': 'gray'}
gl0.ylabel_style = {'size': 12, 'color': 'gray'}

#cax0=axs[0].pcolormesh(lon2d, lat2d, B_p*1.e5, rasterized=True, cmap='magma', transform=trans )
cax0=axs[0].pcolormesh(lon2d, lat2d, B_p*1.e5, vmin=-8., vmax=0., rasterized=True, cmap='magma', transform=trans )

axs[0].pcolormesh(lon2d, lat2d, msk, vmin=-1, vmax=3, transform=trans, rasterized=True, cmap='Greys' )

axs[0].contour( lon2d, lat2d, bathy, [750,1500,2250,3000,3750,4500,5250,6000], 
                colors='gray', alpha=0.7, linewidths=0.6, transform=trans )

axs[0].set_title('(c) Present-day buoyancy flux (1989-2009)',size=14,fontweight='bold')

#colorbar:
#cbar0=fig.colorbar(cax0,ax=axs[0],fraction=0.029, pad=0.02, extend='min')
cbar0=fig.colorbar(cax0,ax=axs[0],fraction=0.029, pad=0.02, extend='min', ticks=np.arange(-8,1,1))
cbar0.ax.set_title('10$^{-5}$ kg.m$^{-1}$.s$^{-3}$',size=12)
cbar0.outline.set_linewidth(0.4)
cbar0.ax.tick_params(labelsize=12,which='both')

# reduce window size:
xmin, xmax = axs[0].get_xlim()
ymin, ymax = axs[0].get_ylim()
axs[0].set_extent( [ xmin + 0.25*(xmax-xmin), \
                     xmax - 0.12*(xmax-xmin), \
                     ymin + 0.14*(ymax-ymin), \
                     ymax ], \
                   crs=proj )

# Fill lower right corner (out of the NEMO domain) ; values from the cursor on plt.show():
x0=8.77e5 ; xlrc=[x0,xmax,xmax,x0]
y0=1.63e6 ; ylrc=[ymin,ymin,y0,y0]
mm=np.zeros((4,4))
axs[0].pcolormesh(xlrc,ylrc,mm,vmin=-1, vmax=3, rasterized=True, cmap='Greys' )

# Plot mean value over the Amundsen continental shelf:
string=np.array2string(mean_B_p.values*1.e5,precision=2,floatmode='fixed')
axs[0].text(xmin+0.365*(xmax-xmin), ymin+0.20*(ymax-ymin),string,color='black',
            fontsize=14,fontweight='bold',ha='right')
axs[0].text(xmin+0.370*(xmax-xmin), ymin+0.20*(ymax-ymin),r'$\times$ 10$^{-5}$ kg.m$^{-1}$.s$^{-3}$',color='black',fontsize=14,fontweight='bold')

#--------------------

#axs[1] = plt.axes(projection=proj)

# (lonmin, lonmax, latmin, latmax) :
axs[1].set_extent( (-145, -85, -76, -68), trans )

# Add and customize meridians/parallels
gl1 = axs[1].gridlines(draw_labels=True, linewidth=0.75, color='gray', alpha=0.7, linestyle='--',
                       dms=True, x_inline=False, y_inline=False )
gl1.top_labels = False
gl1.right_labels = False
gl1.xlocator = mticker.FixedLocator([-140, -130, -120, -110, -100, -90])
gl1.ylocator = mticker.FixedLocator([-78, -76, -74, -72, -70, -68])
gl1.xformatter = LongitudeFormatter()
gl1.yformatter = LatitudeFormatter()
gl1.xlabel_style = {'size': 12, 'color': 'gray'}
gl1.ylabel_style = {'size': 12, 'color': 'gray'}

cax1=axs[1].pcolormesh(lon2d, lat2d, (B_f-B_p)*1.e5, vmin=-3., vmax=3., rasterized=True, cmap='RdBu_r', transform=trans )
#cax1=axs[1].pcolormesh(lon2d, lat2d, (B_f-B_p)*1.e4, rasterized=True, cmap='RdBu_r', transform=trans )

axs[1].pcolormesh(lon2d, lat2d, msk, vmin=-1, vmax=3, transform=trans, rasterized=True, cmap='Greys' )

axs[1].contour( lon2d, lat2d, bathy, [750,1500,2250,3000,3750,4500,5250,6000], 
                colors='gray', alpha=0.7, linewidths=0.6, transform=trans )

axs[1].set_title('(d) Buoyancy flux anomaly',size=14,fontweight='bold')

#colorbar:
#cbar1=fig.colorbar(cax1,ax=axs[1],fraction=0.029, pad=0.02, extend='max')
cbar1=fig.colorbar(cax1,ax=axs[1],fraction=0.029, pad=0.02, extend='both', ticks=np.arange(-3,4,1))
cbar1.ax.set_title('10$^{-5}$ kg.m$^{-1}$.s$^{-3}$ \n ',size=12)
cbar1.outline.set_linewidth(0.4)
cbar1.ax.tick_params(labelsize=12,which='both')

# reduce window size:
xmin, xmax = axs[1].get_xlim()
ymin, ymax = axs[1].get_ylim()
axs[1].set_extent( [ xmin + 0.25*(xmax-xmin), \
                     xmax - 0.12*(xmax-xmin), \
                     ymin + 0.14*(ymax-ymin), \
                     ymax ], \
                   crs=proj )

# Fill lower right corner (out of the NEMO domain) ; values from the cursor on plt.show():
x0=8.77e5 ; xlrc=[x0,xmax,xmax,x0]
y0=1.63e6 ; ylrc=[ymin,ymin,y0,y0]
mm=np.zeros((4,4))
axs[1].pcolormesh(xlrc,ylrc,mm,vmin=-1, vmax=3, rasterized=True, cmap='Greys' )

# Plot mean value over the Amundsen continental shelf:
diff = (mean_B_f - mean_B_p)*1.e5
string=np.array2string(diff.values,precision=2,floatmode='fixed',sign='+')
axs[1].text(xmin+0.365*(xmax-xmin), ymin+0.20*(ymax-ymin),string,color='firebrick',
            fontsize=14,fontweight='bold',ha='right')
axs[1].text(xmin+0.370*(xmax-xmin), ymin+0.20*(ymax-ymin),r'$\times$ 10$^{-5}$ kg.m$^{-1}$.s$^{-3}$',color='firebrick',fontsize=14,fontweight='bold')

#--------------------
#fig.tight_layout()
fig.savefig('map_buoyancy_flux_anomaly.jpg') # warning: do not specify dpi
fig.savefig('map_buoyancy_flux_anomaly.pdf')
#plt.show()
