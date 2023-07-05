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
#msk[msk<0.9]=np.nan
msk[sh<5.0]=np.nan
msk=msk*0.e0+1.e0 # nan over ocean, 1 over ice-sheet (includes ice shelves and grounded ice, but exclude nunatak)
mskb=msk*1.0
mskb[np.isnan(mskb)]=0.0
lon=ncM['LON'].values[j1:j2,i1:i2]
lat=ncM['LAT'].values[j1:j2,i1:i2]
x2d, y2d = np.meshgrid(ncM['x'].values[i1:i2], ncM['y'].values[j1:j2], sparse=False, indexing='xy')

pref_pres='DATA/ERAi_ANx/MARv3.9.3_Amundsen_'
suff_pres='_1979-2017_monthly.nc'
pref_futu='DATA/ERAi_CMIP5anom_ANf/MARv3.9.3_Amundsen_'
suff_futu='_2088-2117_monthly.nc'

file_smb_pres=pref_pres+'smb'+suff_pres
ncSp=xr.open_dataset(file_smb_pres,decode_cf=False)

file_smb_futu=pref_futu+'smb'+suff_futu
ncSf=xr.open_dataset(file_smb_futu,decode_cf=False)

# seasonal means in mm.w.e/(3 months)
smb_anom_SON = np.nansum(   ( ncSf['smb'].values[ 8::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+ 8::12,0,j1:j2,i1:i2] ) * 30. / 30. \
                          + ( ncSf['smb'].values[ 9::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+ 9::12,0,j1:j2,i1:i2] ) * 31. / 30. \
                          + ( ncSf['smb'].values[10::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+10::12,0,j1:j2,i1:i2] ) * 30. / 30. , axis=0) * msk
smb_anom_DJF = np.nansum(   ( ncSf['smb'].values[11::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+11::12,0,j1:j2,i1:i2] ) * 31. / 30. \
                          + ( ncSf['smb'].values[ 0::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+ 0::12,0,j1:j2,i1:i2] ) * 31. / 30. \
                          + ( ncSf['smb'].values[ 1::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+ 1::12,0,j1:j2,i1:i2] ) * 28. / 30. , axis=0) * msk
smb_anom_MAM = np.nansum(   ( ncSf['smb'].values[ 2::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+ 2::12,0,j1:j2,i1:i2] ) * 31. / 30. \
                          + ( ncSf['smb'].values[ 3::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+ 3::12,0,j1:j2,i1:i2] ) * 30. / 30. \
                          + ( ncSf['smb'].values[ 4::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+ 4::12,0,j1:j2,i1:i2] ) * 31. / 30. , axis=0) * msk
smb_anom_JJA = np.nansum(   ( ncSf['smb'].values[ 5::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+ 5::12,0,j1:j2,i1:i2] ) * 30. / 30. \
                          + ( ncSf['smb'].values[ 6::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+ 6::12,0,j1:j2,i1:i2] ) * 31. / 30. \
                          + ( ncSf['smb'].values[ 7::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+ 7::12,0,j1:j2,i1:i2] ) * 31. / 30. , axis=0) * msk

smb_std_SON = np.nanstd(   ( ncSf['smb'].values[ 8::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+ 8::12,0,j1:j2,i1:i2] ) * 30. / 30. \
                         + ( ncSf['smb'].values[ 9::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+ 9::12,0,j1:j2,i1:i2] ) * 31. / 30. \
                         + ( ncSf['smb'].values[10::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+10::12,0,j1:j2,i1:i2] ) * 30. / 30. , axis=0) * msk
smb_std_DJF = np.nanstd(   ( ncSf['smb'].values[11::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+11::12,0,j1:j2,i1:i2] ) * 31. / 30. \
                         + ( ncSf['smb'].values[ 0::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+ 0::12,0,j1:j2,i1:i2] ) * 31. / 30. \
                         + ( ncSf['smb'].values[ 1::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+ 1::12,0,j1:j2,i1:i2] ) * 28. / 30. , axis=0) * msk
smb_std_MAM = np.nanstd(   ( ncSf['smb'].values[ 2::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+ 2::12,0,j1:j2,i1:i2] ) * 31. / 30. \
                         + ( ncSf['smb'].values[ 3::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+ 3::12,0,j1:j2,i1:i2] ) * 30. / 30. \
                         + ( ncSf['smb'].values[ 4::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+ 4::12,0,j1:j2,i1:i2] ) * 31. / 30. , axis=0) * msk
smb_std_JJA = np.nanstd(   ( ncSf['smb'].values[ 5::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+ 5::12,0,j1:j2,i1:i2] ) * 30. / 30. \
                         + ( ncSf['smb'].values[ 6::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+ 6::12,0,j1:j2,i1:i2] ) * 31. / 30. \
                         + ( ncSf['smb'].values[ 7::12,0,j1:j2,i1:i2] - ncSp['smb'].values[9*12+ 7::12,0,j1:j2,i1:i2] ) * 31. / 30. , axis=0) * msk

# t-value
t_SON = smb_anom_SON*np.sqrt(30)/smb_std_SON
t_DJF = smb_anom_DJF*np.sqrt(30)/smb_std_DJF
t_MAM = smb_anom_MAM*np.sqrt(30)/smb_std_MAM
t_JJA = smb_anom_JJA*np.sqrt(30)/smb_std_JJA

# Significant differences (t-test):
t90 = 1.699
t95 = 2.045
t99 = 2.756
tused = t95
print('t value used = ',tused)
eps = 0.1
tmp=t_SON*1.e0
t_SON[np.abs(tmp)>=tused]=1 ; t_SON[np.abs(tmp)<tused]=0 ; t_SON[smb_std_SON<eps]=0
tmp=t_DJF*1.e0
t_DJF[np.abs(tmp)>=tused]=1 ; t_DJF[np.abs(tmp)<tused]=0 ; t_DJF[smb_std_DJF<eps]=0
tmp=t_MAM*1.e0
t_MAM[np.abs(tmp)>=tused]=1 ; t_MAM[np.abs(tmp)<tused]=0 ; t_MAM[smb_std_MAM<eps]=0
tmp=t_JJA*1.e0
t_JJA[np.abs(tmp)>=tused]=1 ; t_JJA[np.abs(tmp)<tused]=0 ; t_JJA[smb_std_JJA<eps]=0

print(np.nanmax(t_SON))
print(np.nanmin(t_SON))
print(np.nanmax(t_JJA))
print(np.nanmin(t_JJA))

#============================================================================
fig, axs = plt.subplots(2,2,figsize=(7.0,7.0))
axs = axs.ravel()
plt.setp(axs, xticks=[], yticks=[])
plt.subplots_adjust(hspace=0.1,bottom=0.2,top=0.8,wspace=0.1)

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

im0=axs[0].contourf(x2d,y2d,smb_anom_SON,np.arange(-150.0,162.5,12.5),cmap=cmap_new,extend='both')
axs[0].contour(x2d,y2d,mskb,[0.5],colors='black',linewidths=0.8)
lo0=axs[0].contour(x2d,y2d,lon,np.arange(-160.,180.,20.),colors='black',linewidths=0.3)
axs[0].clabel(lo0, fmt='%3.0f', colors='black', fontsize=5)
la0=axs[0].contour(x2d,y2d,lat,[-85.,-80.,-75.,-70.],colors='black',linewidths=0.3)
axs[0].clabel(la0, fmt='%2.0f', colors='black', fontsize=5)
sh0=axs[0].contour(x2d,y2d,sh,np.arange(400,4400,400),colors='grey',linewidths=0.5)
axs[0].clabel(sh0, fmt='%4.0f', colors='grey', fontsize=3)
plt.rcParams['hatch.linewidth'] = 0.25
plt.rcParams['hatch.color'] = 'k'
axs[0].contourf(x2d,y2d,t_SON,[-0.5,0.5,1.5],colors='none',hatches=['////////',''])
axs[0].set_aspect(1.0/axs[0].get_data_ratio()*ratio)
axins0 = inset_axes(axs[0], width="90%", height="5%", loc='upper center')
cbar0=fig.colorbar(im0, cax=axins0, orientation="horizontal", ticks=np.arange(-150.0,200.0,50.0))
cbar0.ax.tick_params(labelsize=8)
axins0.xaxis.set_ticks_position("bottom")
axs[0].set_title('(a) SON (mm w. e. (3 months)$^{-1}$)',fontsize=8)

im1=axs[1].contourf(x2d,y2d,smb_anom_DJF,np.arange(-150.0,162.5,12.5),cmap=cmap_new,extend='both')
axs[1].contour(x2d,y2d,mskb,[0.5],colors='black',linewidths=0.8)
lo1=axs[1].contour(x2d,y2d,lon,np.arange(-160.,180.,20.),colors='black',linewidths=0.3)
axs[1].clabel(lo1, fmt='%3.0f', colors='black', fontsize=5)
la1=axs[1].contour(x2d,y2d,lat,[-85.,-80.,-75.,-70.],colors='black',linewidths=0.3)
axs[1].clabel(la1, fmt='%2.0f', colors='black', fontsize=5)
sh1=axs[1].contour(x2d,y2d,sh,np.arange(400,4400,400),colors='grey',linewidths=0.5)
axs[1].clabel(sh1, fmt='%4.0f', colors='grey', fontsize=3)
plt.rcParams['hatch.linewidth'] = 0.25
plt.rcParams['hatch.color'] = 'k'
axs[1].contourf(x2d,y2d,t_DJF,[-0.5,0.5,1.5],colors='none',hatches=['////////',''])
axs[1].set_aspect(1.0/axs[1].get_data_ratio()*ratio)
axins1 = inset_axes(axs[1], width="90%", height="5%", loc='upper center')
cbar1=fig.colorbar(im1, cax=axins1, orientation="horizontal", ticks=np.arange(-150.0,200.0,50.0))
cbar1.ax.tick_params(labelsize=8)
axins1.xaxis.set_ticks_position("bottom")
axs[1].set_title('(b) DJF (mm w. e. (3 months)$^{-1}$)',fontsize=8)

im2=axs[2].contourf(x2d,y2d,smb_anom_MAM,np.arange(-150.0,162.5,12.5),cmap=cmap_new,extend='both')
axs[2].contour(x2d,y2d,mskb,[0.5],colors='black',linewidths=0.8)
lo2=axs[2].contour(x2d,y2d,lon,np.arange(-160.,180.,20.),colors='black',linewidths=0.3)
axs[2].clabel(lo2, fmt='%3.0f', colors='black', fontsize=5)
la2=axs[2].contour(x2d,y2d,lat,[-85.,-80.,-75.,-70.],colors='black',linewidths=0.3)
axs[2].clabel(la2, fmt='%2.0f', colors='black', fontsize=5)
sh2=axs[2].contour(x2d,y2d,sh,np.arange(400,4400,400),colors='grey',linewidths=0.5)
axs[2].clabel(sh2, fmt='%4.0f', colors='grey', fontsize=3)
plt.rcParams['hatch.linewidth'] = 0.25
plt.rcParams['hatch.color'] = 'k'
axs[2].contourf(x2d,y2d,t_MAM,[-0.5,0.5,1.5],colors='none',hatches=['////////',''])
axs[2].set_aspect(1.0/axs[2].get_data_ratio()*ratio)
axins2 = inset_axes(axs[2], width="90%", height="5%", loc='upper center')
cbar2=fig.colorbar(im2, cax=axins2, orientation="horizontal", ticks=np.arange(-150.0,200.0,50.0))
cbar2.ax.tick_params(labelsize=8)
axins2.xaxis.set_ticks_position("bottom")
axs[2].set_title('(c) MAM (mm w. e. (3 months)$^{-1}$)',fontsize=8)

im3=axs[3].contourf(x2d,y2d,smb_anom_JJA,np.arange(-150.0,162.5,12.5),cmap=cmap_new,extend='both')
axs[3].contour(x2d,y2d,mskb,[0.5],colors='black',linewidths=0.8)
lo3=axs[3].contour(x2d,y2d,lon,np.arange(-160.,180.,20.),colors='black',linewidths=0.3)
axs[3].clabel(lo3, fmt='%3.0f', colors='black', fontsize=5)
la3=axs[3].contour(x2d,y2d,lat,[-85.,-80.,-75.,-70.],colors='black',linewidths=0.3)
axs[3].clabel(la3, fmt='%2.0f', colors='black', fontsize=5)
sh3=axs[3].contour(x2d,y2d,sh,np.arange(400,4400,400),colors='grey',linewidths=0.5)
axs[3].clabel(sh3, fmt='%4.0f', colors='grey', fontsize=3)
plt.rcParams['hatch.linewidth'] = 0.25
plt.rcParams['hatch.color'] = 'k'
axs[3].contourf(x2d,y2d,t_JJA,[-0.5,0.5,1.5],colors='none',hatches=['////////',''])
axs[3].set_aspect(1.0/axs[3].get_data_ratio()*ratio)
axins3 = inset_axes(axs[3], width="90%", height="5%", loc='upper center')
cbar3=fig.colorbar(im3, cax=axins3, orientation="horizontal", ticks=np.arange(-150.0,200.0,50.0))
cbar3.ax.tick_params(labelsize=8)
axins3.xaxis.set_ticks_position("bottom")
axs[3].set_title('(d) JJA (mm w. e. (3 months)$^{-1}$)',fontsize=8)

fig.savefig('map_seasonal_SMB_ERAi_CMIP5MMM.pdf')
