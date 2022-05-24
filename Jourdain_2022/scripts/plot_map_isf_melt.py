import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import stereo as st

file_SBC_p_A='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MAR/climato_monthly_AMUXL12-GNJ002_BM02MAR_SBC_1989_2009.nc'
file_SBC_p_B='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM03MAR/climato_monthly_AMUXL12-GNJ002_BM03MAR_SBC_1989_2009.nc'
file_SBC_p_C='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM04MAR/climato_monthly_AMUXL12-GNJ002_BM04MAR_SBC_1989_2009.nc'
file_SBC_f_A='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MARrcp85/climato_monthly_AMUXL12-GNJ002_BM02MARrcp85_SBC_1989_2009.nc'
file_SBC_f_B='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM03MARrcBDY/climato_monthly_AMUXL12-GNJ002_BM03MARrcBDY_SBC_1989_2009.nc'
file_SBC_f_C='/Users/njourdain/OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM04MARrcp85/climato_monthly_AMUXL12-GNJ002_BM04MARrcp85_SBC_1989_2009.nc'

ncpA = xr.open_dataset(file_SBC_p_A,decode_cf=False)
ncpB = xr.open_dataset(file_SBC_p_B,decode_cf=False)
ncpC = xr.open_dataset(file_SBC_p_C,decode_cf=False)
ncfA = xr.open_dataset(file_SBC_f_A,decode_cf=False)
ncfB = xr.open_dataset(file_SBC_f_B,decode_cf=False)
ncfC = xr.open_dataset(file_SBC_f_C,decode_cf=False)

print(ncpA.fwfisf.shape)

# Average over A, B, C simulations (in m.w.e / yr):
isfmelt_p =  - ( ncpA.fwfisf.mean(axis=0) + ncpB.fwfisf.mean(axis=0) + ncpC.fwfisf.mean(axis=0) ) * 1.e-3 * 86400. * 365.25 /3.
isfmelt_f =  - ( ncfA.fwfisf.mean(axis=0) + ncfB.fwfisf.mean(axis=0) + ncfC.fwfisf.mean(axis=0) ) * 1.e-3 * 86400. * 365.25 /3.

file_mshA = '/Users/njourdain/OUTPUT_AMU/mesh_mask_AMUXL12_BedMachineAntarctica-2019-05-24.nc'
file_msh = '/Users/njourdain/OUTPUT_AMU/mesh_mask_AMUXL12_BedMachineAntarctica-2020-07-15_v02.nc'
file_bat = '/Users/njourdain/OUTPUT_AMU/bathy_meter_AMUXL12_BedMachineAntarctica-2020-07-15_v02.nc'

ncmsh = xr.open_dataset(file_msh,decode_cf=False)
ncbat = xr.open_dataset(file_bat,decode_cf=False)
ncmshA = xr.open_dataset(file_mshA,decode_cf=False)

lon2d = ncmsh.glamt.isel(t=0)
lat2d = ncmsh.gphit.isel(t=0)

bathy = ncbat.Bathymetry

msk = ncmsh.tmask.isel(t=0,z=0)*1.0 + ncmsh.tmask.isel(t=0).max('z')*1.0 # 2 open ocean, 1 ice shelf, 0 grounded
msk_grounded = msk.where( (msk<0.5), np.nan ).values
msk_ocean = msk.where( (msk>1.5), np.nan ).values
msk = msk.where( (msk<1.5) & (msk>0.5), np.nan )

mskA = ncmshA.tmask.isel(t=0,z=0)*1.0 + ncmshA.tmask.isel(t=0).max('z')*1.0 # 2 open ocean, 1 ice shelf, 0 grounded
mskA = mskA.where( (mskA<1.5) & (mskA>0.5), np.nan )

#-------------------------------------------------------------------------
# Average over the continental shelf:

msk_isf = ncmsh.tmask.isel(t=0).max('z').where((ncmsh.misf.isel(t=0)>1)&(lon2d>=-135)&(lon2d<=-100.0)&(lat2d<-70.2),0.e0)
tmp_isfmelt_p = isfmelt_p * msk_isf *  ncmsh.e1t.isel(t=0) * ncmsh.e2t.isel(t=0)
tmp_isfmelt_f = isfmelt_f * msk_isf *  ncmsh.e1t.isel(t=0) * ncmsh.e2t.isel(t=0)
total_isfmelt_p = tmp_isfmelt_p.sum() * 1.e-9 # Gt/yr
total_isfmelt_f = tmp_isfmelt_f.sum() * 1.e-9 # Gt/yt

msk_isfA = ncmshA.tmask.isel(t=0).max('z').where((ncmshA.misf.isel(t=0)>1)&(lon2d>=-135)&(lon2d<=-100.0)&(lat2d<-70.2),0.e0)
tmp_isfmelt_p_A = - ncpA.fwfisf.mean(axis=0) * 1.e-3 * 86400. * 365.25 * msk_isfA *  ncmsh.e1t.isel(t=0) * ncmsh.e2t.isel(t=0)
tmp_isfmelt_f_A = - ncfA.fwfisf.mean(axis=0) * 1.e-3 * 86400. * 365.25 * msk_isfA *  ncmsh.e1t.isel(t=0) * ncmsh.e2t.isel(t=0)
total_isfmelt_p_A = tmp_isfmelt_p_A.sum().values * 1.e-9 # Gt/yr
total_isfmelt_f_A = tmp_isfmelt_f_A.sum().values * 1.e-9 # Gt/yt

tmp_isfmelt_p_B = - ncpB.fwfisf.mean(axis=0) * 1.e-3 * 86400. * 365.25 * msk_isf *  ncmsh.e1t.isel(t=0) * ncmsh.e2t.isel(t=0)
tmp_isfmelt_f_B = - ncfB.fwfisf.mean(axis=0) * 1.e-3 * 86400. * 365.25 * msk_isf *  ncmsh.e1t.isel(t=0) * ncmsh.e2t.isel(t=0)
total_isfmelt_p_B = tmp_isfmelt_p_B.sum().values * 1.e-9 # Gt/yr
total_isfmelt_f_B = tmp_isfmelt_f_B.sum().values * 1.e-9 # Gt/yt

tmp_isfmelt_p_C = - ncpC.fwfisf.mean(axis=0) * 1.e-3 * 86400. * 365.25 * msk_isf *  ncmsh.e1t.isel(t=0) * ncmsh.e2t.isel(t=0)
tmp_isfmelt_f_C = - ncfC.fwfisf.mean(axis=0) * 1.e-3 * 86400. * 365.25 * msk_isf *  ncmsh.e1t.isel(t=0) * ncmsh.e2t.isel(t=0)
total_isfmelt_p_C = tmp_isfmelt_p_C.sum().values * 1.e-9 # Gt/yr
total_isfmelt_f_C = tmp_isfmelt_f_C.sum().values * 1.e-9 # Gt/yt

print('total present (A,B,C): ',total_isfmelt_p_A,total_isfmelt_p_B,total_isfmelt_p_C)
print('total future (A,B,C): ',total_isfmelt_f_A,total_isfmelt_f_B,total_isfmelt_f_C)
print('total change (A,B,C): ',total_isfmelt_f_A-total_isfmelt_p_A,total_isfmelt_f_B-total_isfmelt_p_B,total_isfmelt_f_C-total_isfmelt_p_C)
print('mean change: ',np.mean([total_isfmelt_f_A-total_isfmelt_p_A,total_isfmelt_f_B-total_isfmelt_p_B,total_isfmelt_f_C-total_isfmelt_p_C]))
print('std change: ',np.std([total_isfmelt_f_A-total_isfmelt_p_A,total_isfmelt_f_B-total_isfmelt_p_B,total_isfmelt_f_C-total_isfmelt_p_C]))

#-------------------------------------------------------------------------

proj=ccrs.SouthPolarStereo(central_longitude=-115.0)
trans=ccrs.PlateCarree()

fig, axs = plt.subplots(nrows=2,ncols=1,
                        subplot_kw={'projection': proj},
                        figsize=(21.0,12.0))
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
gl0.xlabel_style = {'size': 16, 'color': 'gray'}
gl0.ylabel_style = {'size': 16, 'color': 'gray'}

cax0=axs[0].pcolormesh(lon2d, lat2d, isfmelt_p, vmin=0., vmax=100., rasterized=True, cmap='magma', transform=trans )
#cax0=axs[0].pcolormesh(lon2d, lat2d, isfmelt_p, rasterized=True, cmap='magma', transform=trans )

axs[0].pcolormesh(lon2d, lat2d, msk_grounded, vmin=-1, vmax=3, transform=trans, rasterized=True, cmap='Greys' )
axs[0].pcolormesh(lon2d, lat2d, msk_ocean, vmin=-1, vmax=5, transform=trans, rasterized=True, cmap='Blues' )

axs[0].contour( lon2d, lat2d, bathy, [750,1500,2250,3000,3750,4500,5250,6000], 
                colors='k', alpha=0.7, linewidths=0.6, transform=trans )

axs[0].set_title('(b) Ice shelf basal melt rate (1989-2009)',size=20,fontweight='bold')

#colorbar:
cbar0=fig.colorbar(cax0,ax=axs[0],fraction=0.029, pad=0.02, extend='max', ticks=np.arange(0,110,10))
cbar0.ax.set_title('m.w.e./yr',size=16)
cbar0.outline.set_linewidth(0.4)
cbar0.ax.tick_params(labelsize=16,which='both')

# reduce window size:
xmin, xmax = axs[0].get_xlim()
ymin, ymax = axs[0].get_ylim()
axs[0].set_extent( [ xmin + 0.25*(xmax-xmin),   \
                     xmax - 0.12*(xmax-xmin),   \
                     ymin + 0.18*(ymax-ymin),   \
                     ymax - 0.44*(ymax-ymin) ], \
                   crs=proj )

# Fill lower right corner (out of the NEMO domain) ; values from the cursor on plt.show():
x0=8.77e5 ; xlrc=[x0,xmax,xmax,x0]
y0=1.63e6 ; ylrc=[ymin,ymin,y0,y0]
mm=np.zeros((4,4))
axs[0].pcolormesh(xlrc,ylrc,mm,vmin=-1, vmax=3, rasterized=True, cmap='Greys' )

# Plot mean value over the Amundsen continental shelf:
#string=np.array2string(np.rint(total_isfmelt_p.values))
string=np.char.mod('%d', np.rint(total_isfmelt_p.values))
axs[0].text(xmin+0.365*(xmax-xmin), ymin+0.20*(ymax-ymin),string,color='black',
            fontsize=20,fontweight='bold',ha='right')
axs[0].text(xmin+0.370*(xmax-xmin), ymin+0.20*(ymax-ymin),'Gt/yr',color='black',fontsize=20,fontweight='bold')

axs[0].text(xmin+0.40*(xmax-xmin),ymin+0.260*(ymax-ymin),'Getz',color='black',fontsize=16,fontweight='bold')
axs[0].text(xmin+0.47*(xmax-xmin),ymin+0.275*(ymax-ymin),'Dotson',color='black',fontsize=16,fontweight='bold')
axs[0].text(xmin+0.53*(xmax-xmin),ymin+0.235*(ymax-ymin),'Crosson',color='black',fontsize=16,fontweight='bold')
axs[0].text(xmin+0.585*(xmax-xmin),ymin+0.205*(ymax-ymin),'Thwaites',color='black',fontsize=16,fontweight='bold')
axs[0].text(xmin+0.685*(xmax-xmin),ymin+0.220*(ymax-ymin),'Pine Island',color='black',fontsize=16,fontweight='bold')
axs[0].text(xmin+0.67*(xmax-xmin),ymin+0.345*(ymax-ymin),'Cosgrove',color='black',fontsize=16,fontweight='bold')
axs[0].text(xmin+0.75*(xmax-xmin),ymin+0.320*(ymax-ymin),'Abbot',color='black',fontsize=16,fontweight='bold')
axs[0].text(xmin+0.83*(xmax-xmin),ymin+0.230*(ymax-ymin),'Venable',color='black',fontsize=16,fontweight='bold')



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
gl1.xlabel_style = {'size': 16, 'color': 'gray'}
gl1.ylabel_style = {'size': 16, 'color': 'gray'}

cax1=axs[1].pcolormesh(lon2d, lat2d, isfmelt_f, vmin=0., vmax=100., rasterized=True, cmap='magma', transform=trans )
#cax1=axs[1].pcolormesh(lon2d, lat2d, isfmelt_f, rasterized=True, cmap='magma', transform=trans )

axs[1].pcolormesh(lon2d, lat2d, msk_grounded, vmin=-1, vmax=3, transform=trans, rasterized=True, cmap='Greys' )
axs[1].pcolormesh(lon2d, lat2d, msk_ocean, vmin=-1, vmax=5, transform=trans, rasterized=True, cmap='Blues' )

axs[1].contour( lon2d, lat2d, bathy, [750,1500,2250,3000,3750,4500,5250,6000], 
                colors='k', alpha=0.7, linewidths=0.6, transform=trans )

axs[1].set_title('(c) Ice shelf basal melt rate (2080-2100)',size=20,fontweight='bold')

#colorbar:
cbar1=fig.colorbar(cax1,ax=axs[1], fraction=0.029, pad=0.02, extend='max', ticks=np.arange(0,110,10))
cbar1.ax.set_title('m.w.e./yr',size=16)
cbar1.outline.set_linewidth(0.4)
cbar1.ax.tick_params(labelsize=16,which='both')

# reduce window size:
xmin, xmax = axs[1].get_xlim()
ymin, ymax = axs[1].get_ylim()
axs[1].set_extent( [ xmin + 0.25*(xmax-xmin),   \
                     xmax - 0.12*(xmax-xmin),   \
                     ymin + 0.18*(ymax-ymin),   \
                     ymax - 0.44*(ymax-ymin) ], \
                   crs=proj )

# Fill lower right corner (out of the NEMO domain) ; values from the cursor on plt.show():
x0=8.77e5 ; xlrc=[x0,xmax,xmax,x0]
y0=1.63e6 ; ylrc=[ymin,ymin,y0,y0]
mm=np.zeros((4,4))
axs[1].pcolormesh(xlrc,ylrc,mm,vmin=-1, vmax=3, rasterized=True, cmap='Greys' )

# Plot mean value over the Amundsen continental shelf:
string=np.char.mod('%d', np.rint(total_isfmelt_f.values))
axs[1].text(xmin+0.365*(xmax-xmin), ymin+0.20*(ymax-ymin),string,color='k',
            fontsize=20,fontweight='bold',ha='right')
axs[1].text(xmin+0.370*(xmax-xmin), ymin+0.20*(ymax-ymin),'Gt/yr',color='k',fontsize=20,fontweight='bold')


#--------------------
fig.savefig('map_isf_melt.jpg') # warning: do not specify dpi
fig.savefig('map_isf_melt.pdf')
#plt.show()
