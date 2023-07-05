import numpy as np
import xarray as xr
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Domain reduction for plots:
i1=6
i2=-6
j1=6
j2=-45

file_msk='MAR_grid10km.nc'
ncM=xr.open_dataset(file_msk,decode_cf=False)
sh=ncM['SH'].values[j1:j2,i1:i2]
msk=ncM['MSK'].values[j1:j2,i1:i2]
msk[sh<5.0]=np.nan
msk=msk*0.e0+1.e0 # nan over ocean, 1 over ice-sheet (includes ice shelves and grounded ice, but exclude nunatak)
mskb=msk*1.0
mskb[np.isnan(mskb)]=0.0
lon=ncM['LON'].values[j1:j2,i1:i2]
lat=ncM['LAT'].values[j1:j2,i1:i2]
x2d, y2d = np.meshgrid(ncM['x'].values[i1:i2], ncM['y'].values[j1:j2], sparse=False, indexing='xy')

msk_isf=sh*1
msk_isf[sh<10.0]=np.nan
msk_isf[sh>200.0]=np.nan
msk_isf[msk_isf>1]=1
msk_isf[np.isnan(msk_isf)]=0

msk_ABB=np.load('ANT_MASK_ABBOT.npy');   msk_ABB[ np.isnan(msk_ABB) ]=0 ; msk_ABB=msk_ABB[j1:j2,i1:i2]
msk_COS=np.load('ANT_MASK_COSGR.npy');   msk_COS[ np.isnan(msk_COS) ]=0 ; msk_COS=msk_COS[j1:j2,i1:i2]
msk_PIG=np.load('ANT_MASK_PINE.npy');    msk_PIG[ np.isnan(msk_PIG) ]=0 ; msk_PIG=msk_PIG[j1:j2,i1:i2]
msk_THW=np.load('ANT_MASK_THWAIT.npy');  msk_THW[ np.isnan(msk_THW) ]=0 ; msk_THW=msk_THW[j1:j2,i1:i2]
msk_CRO=np.load('ANT_MASK_CROSSON.npy'); msk_CRO[ np.isnan(msk_CRO) ]=0 ; msk_CRO=msk_CRO[j1:j2,i1:i2]
msk_DOT=np.load('ANT_MASK_DOTSON.npy');  msk_DOT[ np.isnan(msk_DOT) ]=0 ; msk_DOT=msk_DOT[j1:j2,i1:i2]
msk_GET=np.load('ANT_MASK_GETZ.npy');    msk_GET[ np.isnan(msk_GET) ]=0 ; msk_GET=msk_GET[j1:j2,i1:i2]

pref_pres='DATA/ACCESS_pres_ANg/MARv3.9.3_Amundsen_'
suff_pres='_1989-2009_ACCESS_monthly.nc'
pref_futu='DATA/ACCESS_futu_ANh/MARv3.9.3_Amundsen_'
suff_futu='_2080-2100_ACCESS_monthly.nc'
pref_anom='DATA/ACCESS_futu_anom_withSICcor_ANm/MARv3.9.3_Amundsen_'
suff_anom='_1989-2009_ANOM_ACCESS_monthly.nc'

file_smb=pref_pres+'smb'+suff_pres
ncSp=xr.open_dataset(file_smb,decode_cf=False)
smb_ACCESS_pres=np.nanmean(ncSp['smb'].values[:,0,j1:j2,i1:i2],axis=0)*msk*365.25 # in mm.w.e/yr

file_smb=pref_futu+'smb'+suff_futu
ncSf=xr.open_dataset(file_smb,decode_cf=False)
smb_ACCESS_futu=np.nanmean(ncSf['smb'].values[:,0,j1:j2,i1:i2],axis=0)*msk*365.25

file_smb=pref_anom+'smb'+suff_anom
ncSa=xr.open_dataset(file_smb,decode_cf=False)
smb_ACCESS_anom=np.nanmean(ncSa['smb'].values[:,0,j1:j2,i1:i2],axis=0)*msk*365.25

file_mlt=pref_pres+'mlt'+suff_pres
ncSp=xr.open_dataset(file_mlt,decode_cf=False)
mlt_ACCESS_pres=np.nanmean(ncSp['mlt'].values[:,0,j1:j2,i1:i2],axis=0)*msk*365.25 # in mm.w.e/yr

file_mlt=pref_futu+'mlt'+suff_futu
ncSf=xr.open_dataset(file_mlt,decode_cf=False)
mlt_ACCESS_futu=np.nanmean(ncSf['mlt'].values[:,0,j1:j2,i1:i2],axis=0)*msk*365.25

file_mlt=pref_anom+'mlt'+suff_anom
ncSa=xr.open_dataset(file_mlt,decode_cf=False)
mlt_ACCESS_anom=np.nanmean(ncSa['mlt'].values[:,0,j1:j2,i1:i2],axis=0)*msk*365.25

#============================================================================
fig, axs = plt.subplots(2,2,figsize=(7.0,7.0))
axs = axs.ravel()
plt.setp(axs, xticks=[], yticks=[])
plt.subplots_adjust(hspace=0.1,bottom=0.2,top=0.8,wspace=0.15)

#----------
# Defining colormap:

# moving the zero of colorbar
# NB: modify the Ncool to Nwarm ratio (total=256) to place zero as desired.
Ncool=128
Nwarm=256-Ncool
col = cm.get_cmap('PuOr', 256)
tmp1 = col(np.linspace(0.47, 1.00, Ncool)) # decrease first number to have more white in the middle light-blue colors
tmp2 = col(np.linspace(0.00, 0.51, Nwarm)) # increase second number to have more white in the middle light-yellow colors
newcolors = np.append(tmp1[::-1,:],tmp2[::-1,:],axis=0)
cmap_new = ListedColormap(newcolors)

#----------
ratio=np.shape(sh)[0]*1./(1.*np.shape(sh)[1])
print(ratio)

im0=axs[0].contourf(x2d,y2d,smb_ACCESS_futu-smb_ACCESS_pres,np.arange(-400.0,440.0,40.0),cmap=cmap_new,extend='both')
axs[0].contour(x2d,y2d,mskb,[0.5],colors='black',linewidths=0.8)
lo0=axs[0].contour(x2d,y2d,lon,np.arange(-160.,180.,20.),colors='black',linewidths=0.3)
axs[0].clabel(lo0, fmt='%3.0f', colors='black', fontsize=5)
la0=axs[0].contour(x2d,y2d,lat,[-85.,-80.,-75.,-70.],colors='black',linewidths=0.3)
axs[0].clabel(la0, fmt='%2.0f', colors='black', fontsize=5)
sh0=axs[0].contour(x2d,y2d,sh,np.arange(400,4400,400),colors='grey',linewidths=0.5)
axs[0].clabel(sh0, fmt='%4.0f', colors='grey', fontsize=3)
axs[0].contour(x2d,y2d,msk_ABB,[0.5],colors='firebrick',linewidths=0.9)
axs[0].text(600.0,220.0,'ABB',weight='bold',color='firebrick',fontsize=5)
axs[0].contour(x2d,y2d,msk_COS,[0.5],colors='firebrick',linewidths=0.9)
axs[0].text(40.0,170.0,'COS',weight='bold',color='firebrick',fontsize=5)
axs[0].plot([210.,350.],[170.,130.],color='firebrick',linewidth=0.7)
axs[0].contour(x2d,y2d,msk_PIG,[0.5],colors='firebrick',linewidths=0.9)
axs[0].text(450.0,-100.0,'PIG',weight='bold',color='firebrick',fontsize=5)
axs[0].contour(x2d,y2d,msk_THW,[0.5],colors='firebrick',linewidths=0.9)
axs[0].text(0.0,-350.0,'THW',weight='bold',color='firebrick',fontsize=5)
axs[0].contour(x2d,y2d,msk_CRO,[0.5],colors='firebrick',linewidths=0.9)
axs[0].text(60.0,80.0,'CRO',weight='bold',color='firebrick',fontsize=5)
axs[0].plot([120.,20.],[60.,-30.],color='firebrick',linewidth=0.7)
axs[0].contour(x2d,y2d,msk_DOT,[0.5],colors='firebrick',linewidths=0.9)
axs[0].text(-170.0,250.0,'DOT',weight='bold',color='firebrick',fontsize=5)
axs[0].plot([-70.,-20.],[220.,50.],color='firebrick',linewidth=0.7)
axs[0].contour(x2d,y2d,msk_GET,[0.5],colors='firebrick',linewidths=0.9)
axs[0].text(-700.0,0.0,'GET',weight='bold',color='firebrick',fontsize=5)
axs[0].set_aspect(1.0/axs[0].get_data_ratio()*ratio)
axins0 = inset_axes(axs[0], width="90%", height="5%", loc='upper center')
cbar0=fig.colorbar(im0, cax=axins0, orientation="horizontal", ticks=np.arange(-400.0,600.0,200.0))
cbar0.ax.tick_params(labelsize=8)
axins0.xaxis.set_ticks_position("bottom")
axs[0].set_title('(a) Future (ACCESS) - Present (ACCESS) \n SMB (mm w. e. yr$^{-1}$)',fontsize=8)

im1=axs[1].contourf(x2d,y2d,smb_ACCESS_futu-smb_ACCESS_anom,np.arange(-400.0,440.0,40.0),cmap=cmap_new,extend='both')
axs[1].contour(x2d,y2d,mskb,[0.5],colors='black',linewidths=0.8)
lo1=axs[1].contour(x2d,y2d,lon,np.arange(-160.,180.,20.),colors='black',linewidths=0.3)
axs[1].clabel(lo1, fmt='%3.0f', colors='black', fontsize=5)
la1=axs[1].contour(x2d,y2d,lat,[-85.,-80.,-75.,-70.],colors='black',linewidths=0.3)
axs[1].clabel(la1, fmt='%2.0f', colors='black', fontsize=5)
sh1=axs[1].contour(x2d,y2d,sh,np.arange(400,4400,400),colors='grey',linewidths=0.5)
axs[1].clabel(sh1, fmt='%4.0f', colors='grey', fontsize=3)
axs[1].set_aspect(1.0/axs[1].get_data_ratio()*ratio)
axins1 = inset_axes(axs[1], width="90%", height="5%", loc='upper center')
cbar1=fig.colorbar(im1, cax=axins1, orientation="horizontal", ticks=np.arange(-400.0,600.0,200.0))
cbar1.ax.tick_params(labelsize=8)
axins1.xaxis.set_ticks_position("bottom")
axs[1].set_title('(b) Future (ACCESS) - Future (anomaly) \n SMB (mm w. e. yr$^{-1}$)',fontsize=8)

im2=axs[2].contourf(x2d,y2d,mlt_ACCESS_futu-mlt_ACCESS_pres,np.arange(-750.0,825.0,75.0),cmap=cmap_new,extend='both')
axs[2].contour(x2d,y2d,mskb,[0.5],colors='black',linewidths=0.8)
lo2=axs[2].contour(x2d,y2d,lon,np.arange(-160.,180.,20.),colors='black',linewidths=0.3)
axs[2].clabel(lo2, fmt='%3.0f', colors='black', fontsize=5)
la2=axs[2].contour(x2d,y2d,lat,[-85.,-80.,-75.,-70.],colors='black',linewidths=0.3)
axs[2].clabel(la2, fmt='%2.0f', colors='black', fontsize=5)
sh2=axs[2].contour(x2d,y2d,sh,np.arange(400,4400,400),colors='grey',linewidths=0.5)
axs[2].clabel(sh2, fmt='%4.0f', colors='grey', fontsize=3)
axs[2].text(650.0,270.0,'ABB',weight='bold',color='firebrick',fontsize=5)
axs[2].plot([550.,640.],[190.,270.],color='firebrick',linewidth=0.7)
axs[2].text(40.0,170.0,'COS',weight='bold',color='firebrick',fontsize=5)
axs[2].plot([210.,350.],[170.,130.],color='firebrick',linewidth=0.7)
axs[2].text(450.0,-120.0,'PIG',weight='bold',color='firebrick',fontsize=5)
axs[2].plot([430.,320.],[-80.,-20.],color='firebrick',linewidth=0.7)
axs[2].text(0.0,-350.0,'THW',weight='bold',color='firebrick',fontsize=5)
axs[2].plot([100.,150.],[-270.,-20.],color='firebrick',linewidth=0.7)
axs[2].text(60.0,80.0,'CRO',weight='bold',color='firebrick',fontsize=5)
axs[2].plot([120.,60.],[70.,0.],color='firebrick',linewidth=0.7)
axs[2].text(-170.0,250.0,'DOT',weight='bold',color='firebrick',fontsize=5)
axs[2].plot([-70.,-20.],[220.,50.],color='firebrick',linewidth=0.7)
axs[2].text(-730.0,0.0,'GET',weight='bold',color='firebrick',fontsize=5)
axs[2].plot([-570.,-420.],[0.,-30.],color='firebrick',linewidth=0.7)
plt.rcParams['hatch.linewidth'] = 0.25
plt.rcParams['hatch.color'] = 'firebrick'
axs[2].contourf(x2d,y2d,msk_isf,[0.0,0.5,1.0],hatches=['','\\\\\\\\\\'],colors='none')
axs[2].set_aspect(1.0/axs[2].get_data_ratio()*ratio)
axins2 = inset_axes(axs[2], width="90%", height="5%", loc='upper center')
cbar2=fig.colorbar(im2, cax=axins2, orientation="horizontal", ticks=np.arange(-600.0,900.0,300.0))
cbar2.ax.tick_params(labelsize=8)
axins2.xaxis.set_ticks_position("bottom")
axs[2].set_title('(c) Future (ACCESS) - Present (ACCESS) \n melt rate (mm w. e. yr$^{-1}$)',fontsize=8)

im3=axs[3].contourf(x2d,y2d,mlt_ACCESS_futu-mlt_ACCESS_anom,np.arange(-750.0,825.0,75.0),cmap=cmap_new,extend='both')
axs[3].contour(x2d,y2d,mskb,[0.5],colors='black',linewidths=0.8)
lo3=axs[3].contour(x2d,y2d,lon,np.arange(-160.,180.,20.),colors='black',linewidths=0.3)
axs[3].clabel(lo3, fmt='%3.0f', colors='black', fontsize=5)
la3=axs[3].contour(x2d,y2d,lat,[-85.,-80.,-75.,-70.],colors='black',linewidths=0.3)
axs[3].clabel(la3, fmt='%2.0f', colors='black', fontsize=5)
sh3=axs[3].contour(x2d,y2d,sh,np.arange(400,4400,400),colors='grey',linewidths=0.5)
axs[3].clabel(sh3, fmt='%4.0f', colors='grey', fontsize=3)
axs[3].set_aspect(1.0/axs[3].get_data_ratio()*ratio)
axins3 = inset_axes(axs[3], width="90%", height="5%", loc='upper center')
cbar3=fig.colorbar(im3, cax=axins3, orientation="horizontal", ticks=np.arange(-600.0,900.0,300.0))
cbar3.ax.tick_params(labelsize=8)
axins3.xaxis.set_ticks_position("bottom")
axs[3].set_title('(d) Future (ACCESS) - Future (anomaly) \n melt rate (mm w. e. yr$^{-1}$)',fontsize=8)

fig.savefig('map_mean_SMB_MLT_ACCESS.pdf')
