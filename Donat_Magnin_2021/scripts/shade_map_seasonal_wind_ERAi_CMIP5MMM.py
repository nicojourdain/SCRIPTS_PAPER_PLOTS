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

file_UUz_pres=pref_pres+'UUz'+suff_pres
ncUp=xr.open_dataset(file_UUz_pres,decode_cf=False)
aa=ncUp['zuvlev'].values
k10=np.unravel_index(np.argmin(np.abs(aa-10.0), axis=None), aa.shape)[0]
print(k10)

file_UUz_futu=pref_futu+'UUz'+suff_futu
ncUf=xr.open_dataset(file_UUz_futu,decode_cf=False)

# seasonal means in m/s
UUz_anom_SON = np.nanmean(   ( ncUf['UUz'].values[ 8::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+ 8::12,k10,j1:j2,i1:i2] ) * 30. / 30. \
                           + ( ncUf['UUz'].values[ 9::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+ 9::12,k10,j1:j2,i1:i2] ) * 31. / 30. \
                           + ( ncUf['UUz'].values[10::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+10::12,k10,j1:j2,i1:i2] ) * 30. / 30. , axis=0) / 3.e0
UUz_anom_DJF = np.nanmean(   ( ncUf['UUz'].values[11::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+11::12,k10,j1:j2,i1:i2] ) * 31. / 30. \
                           + ( ncUf['UUz'].values[ 0::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+ 0::12,k10,j1:j2,i1:i2] ) * 31. / 30. \
                           + ( ncUf['UUz'].values[ 1::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+ 1::12,k10,j1:j2,i1:i2] ) * 28. / 30. , axis=0) / 3.e0
UUz_anom_MAM = np.nanmean(   ( ncUf['UUz'].values[ 2::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+ 2::12,k10,j1:j2,i1:i2] ) * 31. / 30. \
                           + ( ncUf['UUz'].values[ 3::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+ 3::12,k10,j1:j2,i1:i2] ) * 30. / 30. \
                           + ( ncUf['UUz'].values[ 4::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+ 4::12,k10,j1:j2,i1:i2] ) * 31. / 30. , axis=0) / 3.e0
UUz_anom_JJA = np.nanmean(   ( ncUf['UUz'].values[ 5::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+ 5::12,k10,j1:j2,i1:i2] ) * 30. / 30. \
                           + ( ncUf['UUz'].values[ 6::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+ 6::12,k10,j1:j2,i1:i2] ) * 31. / 30. \
                           + ( ncUf['UUz'].values[ 7::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+ 7::12,k10,j1:j2,i1:i2] ) * 31. / 30. , axis=0) / 3.e0

UUz_std_SON = np.nanstd(   ( ncUf['UUz'].values[ 8::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+ 8::12,k10,j1:j2,i1:i2] ) * 30. / 30. \
                         + ( ncUf['UUz'].values[ 9::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+ 9::12,k10,j1:j2,i1:i2] ) * 31. / 30. \
                         + ( ncUf['UUz'].values[10::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+10::12,k10,j1:j2,i1:i2] ) * 30. / 30. , axis=0) / 3.e0
UUz_std_DJF = np.nanstd(   ( ncUf['UUz'].values[11::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+11::12,k10,j1:j2,i1:i2] ) * 31. / 30. \
                         + ( ncUf['UUz'].values[ 0::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+ 0::12,k10,j1:j2,i1:i2] ) * 31. / 30. \
                         + ( ncUf['UUz'].values[ 1::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+ 1::12,k10,j1:j2,i1:i2] ) * 28. / 30. , axis=0) / 3.e0
UUz_std_MAM = np.nanstd(   ( ncUf['UUz'].values[ 2::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+ 2::12,k10,j1:j2,i1:i2] ) * 31. / 30. \
                         + ( ncUf['UUz'].values[ 3::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+ 3::12,k10,j1:j2,i1:i2] ) * 30. / 30. \
                         + ( ncUf['UUz'].values[ 4::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+ 4::12,k10,j1:j2,i1:i2] ) * 31. / 30. , axis=0) / 3.e0
UUz_std_JJA = np.nanstd(   ( ncUf['UUz'].values[ 5::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+ 5::12,k10,j1:j2,i1:i2] ) * 30. / 30. \
                         + ( ncUf['UUz'].values[ 6::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+ 6::12,k10,j1:j2,i1:i2] ) * 31. / 30. \
                         + ( ncUf['UUz'].values[ 7::12,k10,j1:j2,i1:i2] - ncUp['UUz'].values[9*12+ 7::12,k10,j1:j2,i1:i2] ) * 31. / 30. , axis=0) / 3.e0

# t-value
tUU_SON = UUz_anom_SON*np.sqrt(30)/UUz_std_SON
tUU_DJF = UUz_anom_DJF*np.sqrt(30)/UUz_std_DJF
tUU_MAM = UUz_anom_MAM*np.sqrt(30)/UUz_std_MAM
tUU_JJA = UUz_anom_JJA*np.sqrt(30)/UUz_std_JJA

# Significant differences (t-test):
t90 = 1.699
t95 = 2.045
t99 = 2.756
tused = t95
print('t value used = ',tused)
eps = 0.001
tmp=tUU_SON*1.e0
tUU_SON[np.abs(tmp)>=tused]=1 ; tUU_SON[np.abs(tmp)<tused]=0 ; tUU_SON[UUz_std_SON<eps]=0
tmp=tUU_DJF*1.e0
tUU_DJF[np.abs(tmp)>=tused]=1 ; tUU_DJF[np.abs(tmp)<tused]=0 ; tUU_DJF[UUz_std_DJF<eps]=0
tmp=tUU_MAM*1.e0
tUU_MAM[np.abs(tmp)>=tused]=1 ; tUU_MAM[np.abs(tmp)<tused]=0 ; tUU_MAM[UUz_std_MAM<eps]=0
tmp=tUU_JJA*1.e0
tUU_JJA[np.abs(tmp)>=tused]=1 ; tUU_JJA[np.abs(tmp)<tused]=0 ; tUU_JJA[UUz_std_JJA<eps]=0

print(np.nanmax(tUU_SON))
print(np.nanmin(tUU_SON))
print(np.nanmax(tUU_JJA))
print(np.nanmin(tUU_JJA))

#-----

file_VVz_pres=pref_pres+'VVz'+suff_pres
ncVp=xr.open_dataset(file_VVz_pres,decode_cf=False)

file_VVz_futu=pref_futu+'VVz'+suff_futu
ncVf=xr.open_dataset(file_VVz_futu,decode_cf=False)

# seasonal means in m/s
VVz_anom_SON = np.nanmean(   ( ncVf['VVz'].values[ 8::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+ 8::12,k10,j1:j2,i1:i2] ) * 30. / 30. \
                           + ( ncVf['VVz'].values[ 9::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+ 9::12,k10,j1:j2,i1:i2] ) * 31. / 30. \
                           + ( ncVf['VVz'].values[10::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+10::12,k10,j1:j2,i1:i2] ) * 30. / 30. , axis=0) / 3.e0
VVz_anom_DJF = np.nanmean(   ( ncVf['VVz'].values[11::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+11::12,k10,j1:j2,i1:i2] ) * 31. / 30. \
                           + ( ncVf['VVz'].values[ 0::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+ 0::12,k10,j1:j2,i1:i2] ) * 31. / 30. \
                           + ( ncVf['VVz'].values[ 1::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+ 1::12,k10,j1:j2,i1:i2] ) * 28. / 30. , axis=0) / 3.e0
VVz_anom_MAM = np.nanmean(   ( ncVf['VVz'].values[ 2::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+ 2::12,k10,j1:j2,i1:i2] ) * 31. / 30. \
                           + ( ncVf['VVz'].values[ 3::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+ 3::12,k10,j1:j2,i1:i2] ) * 30. / 30. \
                           + ( ncVf['VVz'].values[ 4::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+ 4::12,k10,j1:j2,i1:i2] ) * 31. / 30. , axis=0) / 3.e0
VVz_anom_JJA = np.nanmean(   ( ncVf['VVz'].values[ 5::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+ 5::12,k10,j1:j2,i1:i2] ) * 30. / 30. \
                           + ( ncVf['VVz'].values[ 6::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+ 6::12,k10,j1:j2,i1:i2] ) * 31. / 30. \
                           + ( ncVf['VVz'].values[ 7::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+ 7::12,k10,j1:j2,i1:i2] ) * 31. / 30. , axis=0) / 3.e0

VVz_std_SON = np.nanstd(   ( ncVf['VVz'].values[ 8::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+ 8::12,k10,j1:j2,i1:i2] ) * 30. / 30. \
                         + ( ncVf['VVz'].values[ 9::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+ 9::12,k10,j1:j2,i1:i2] ) * 31. / 30. \
                         + ( ncVf['VVz'].values[10::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+10::12,k10,j1:j2,i1:i2] ) * 30. / 30. , axis=0) / 3.e0
VVz_std_DJF = np.nanstd(   ( ncVf['VVz'].values[11::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+11::12,k10,j1:j2,i1:i2] ) * 31. / 30. \
                         + ( ncVf['VVz'].values[ 0::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+ 0::12,k10,j1:j2,i1:i2] ) * 31. / 30. \
                         + ( ncVf['VVz'].values[ 1::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+ 1::12,k10,j1:j2,i1:i2] ) * 28. / 30. , axis=0) / 3.e0
VVz_std_MAM = np.nanstd(   ( ncVf['VVz'].values[ 2::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+ 2::12,k10,j1:j2,i1:i2] ) * 31. / 30. \
                         + ( ncVf['VVz'].values[ 3::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+ 3::12,k10,j1:j2,i1:i2] ) * 30. / 30. \
                         + ( ncVf['VVz'].values[ 4::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+ 4::12,k10,j1:j2,i1:i2] ) * 31. / 30. , axis=0) / 3.e0
VVz_std_JJA = np.nanstd(   ( ncVf['VVz'].values[ 5::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+ 5::12,k10,j1:j2,i1:i2] ) * 30. / 30. \
                         + ( ncVf['VVz'].values[ 6::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+ 6::12,k10,j1:j2,i1:i2] ) * 31. / 30. \
                         + ( ncVf['VVz'].values[ 7::12,k10,j1:j2,i1:i2] - ncVp['VVz'].values[9*12+ 7::12,k10,j1:j2,i1:i2] ) * 31. / 30. , axis=0) / 3.e0

# t-value
tVV_SON = VVz_anom_SON*np.sqrt(30)/VVz_std_SON
tVV_DJF = VVz_anom_DJF*np.sqrt(30)/VVz_std_DJF
tVV_MAM = VVz_anom_MAM*np.sqrt(30)/VVz_std_MAM
tVV_JJA = VVz_anom_JJA*np.sqrt(30)/VVz_std_JJA

# Significant differences (t-test):
tmp=tVV_SON*1.e0
tVV_SON[np.abs(tmp)>=tused]=1 ; tVV_SON[np.abs(tmp)<tused]=0 ; tVV_SON[VVz_std_SON<eps]=0
tmp=tVV_DJF*1.e0
tVV_DJF[np.abs(tmp)>=tused]=1 ; tVV_DJF[np.abs(tmp)<tused]=0 ; tVV_DJF[VVz_std_DJF<eps]=0
tmp=tVV_MAM*1.e0
tVV_MAM[np.abs(tmp)>=tused]=1 ; tVV_MAM[np.abs(tmp)<tused]=0 ; tVV_MAM[VVz_std_MAM<eps]=0
tmp=tVV_JJA*1.e0
tVV_JJA[np.abs(tmp)>=tused]=1 ; tVV_JJA[np.abs(tmp)<tused]=0 ; tVV_JJA[VVz_std_JJA<eps]=0

print(np.nanmax(tVV_SON))
print(np.nanmin(tVV_SON))
print(np.nanmax(tVV_JJA))
print(np.nanmin(tVV_JJA))

UUz_anom_SON[((tUU_SON==0)|(tVV_SON==0))]=np.nan
UUz_anom_DJF[((tUU_DJF==0)|(tVV_DJF==0))]=np.nan
UUz_anom_MAM[((tUU_MAM==0)|(tVV_MAM==0))]=np.nan
UUz_anom_JJA[((tUU_JJA==0)|(tVV_JJA==0))]=np.nan
VVz_anom_SON[((tUU_SON==0)|(tVV_SON==0))]=np.nan
VVz_anom_DJF[((tUU_DJF==0)|(tVV_DJF==0))]=np.nan
VVz_anom_MAM[((tUU_MAM==0)|(tVV_MAM==0))]=np.nan
VVz_anom_JJA[((tUU_JJA==0)|(tVV_JJA==0))]=np.nan

#============================================================================
fig, axs = plt.subplots(2,2,figsize=(7.0,7.0))
axs = axs.ravel()
plt.setp(axs, xticks=[], yticks=[])
plt.subplots_adjust(hspace=0.1,bottom=0.2,top=0.8,wspace=0.1)

#----------
ratio=np.shape(sh)[0]*1./(1.*np.shape(sh)[1])
print(ratio)

im0=axs[0].contourf(x2d,y2d,mskb,[-0.2,-0.1,0.,0.1,0.2,0.3],cmap='Blues')
axs[0].contour(x2d,y2d,mskb,[0.5],colors='black',linewidths=0.8)
lo0=axs[0].contour(x2d,y2d,lon,np.arange(-160.,180.,20.),colors='black',linewidths=0.3)
la0=axs[0].contour(x2d,y2d,lat,[-85.,-80.,-75.,-70.],colors='black',linewidths=0.3)
sh0=axs[0].contour(x2d,y2d,sh,np.arange(400,4400,400),colors='grey',linewidths=0.5)
q0=axs[0].quiver(x2d[::5,::5],y2d[::5,::5],UUz_anom_SON[::5,::5],VVz_anom_SON[::5,::5],units='width',scale=1/0.05,scale_units='width')
axs[0].quiverkey(q0, 0.78, 0.19, 1.0, '1.0 m s$^{-1}$', labelpos='E',coordinates='figure',fontproperties={'size':7})
axs[0].set_aspect(1.0/axs[0].get_data_ratio()*ratio)
axs[0].set_title('(a) SON',fontsize=8)

im1=axs[1].contourf(x2d,y2d,mskb,[-0.2,-0.1,0.,0.1,0.2,0.3],cmap='Blues')
axs[1].contour(x2d,y2d,mskb,[0.5],colors='black',linewidths=0.8)
lo1=axs[1].contour(x2d,y2d,lon,np.arange(-160.,180.,20.),colors='black',linewidths=0.3)
la1=axs[1].contour(x2d,y2d,lat,[-85.,-80.,-75.,-70.],colors='black',linewidths=0.3)
sh1=axs[1].contour(x2d,y2d,sh,np.arange(400,4400,400),colors='grey',linewidths=0.5)
q1=axs[1].quiver(x2d[::5,::5],y2d[::5,::5],UUz_anom_DJF[::5,::5],VVz_anom_DJF[::5,::5],units='width',scale=1/0.05,scale_units='width')
axs[1].set_aspect(1.0/axs[1].get_data_ratio()*ratio)
axs[1].set_title('(b) DJF',fontsize=8)

im2=axs[2].contourf(x2d,y2d,mskb,[-0.2,-0.1,0.,0.1,0.2,0.3],cmap='Blues')
axs[2].contour(x2d,y2d,mskb,[0.5],colors='black',linewidths=0.8)
lo2=axs[2].contour(x2d,y2d,lon,np.arange(-160.,180.,20.),colors='black',linewidths=0.3)
la2=axs[2].contour(x2d,y2d,lat,[-85.,-80.,-75.,-70.],colors='black',linewidths=0.3)
sh2=axs[2].contour(x2d,y2d,sh,np.arange(400,4400,400),colors='grey',linewidths=0.5)
q2=axs[2].quiver(x2d[::5,::5],y2d[::5,::5],UUz_anom_MAM[::5,::5],VVz_anom_MAM[::5,::5],units='width',scale=1/0.05,scale_units='width')
axs[2].set_aspect(1.0/axs[2].get_data_ratio()*ratio)
axs[2].set_title('(c) MAM',fontsize=8)

im3=axs[3].contourf(x2d,y2d,mskb,[-0.2,-0.1,0.,0.1,0.2,0.3],cmap='Blues')
axs[3].contour(x2d,y2d,mskb,[0.5],colors='black',linewidths=0.8)
lo3=axs[3].contour(x2d,y2d,lon,np.arange(-160.,180.,20.),colors='black',linewidths=0.3)
la3=axs[3].contour(x2d,y2d,lat,[-85.,-80.,-75.,-70.],colors='black',linewidths=0.3)
sh3=axs[3].contour(x2d,y2d,sh,np.arange(400,4400,400),colors='grey',linewidths=0.5)
q3=axs[3].quiver(x2d[::5,::5],y2d[::5,::5],UUz_anom_JJA[::5,::5],VVz_anom_JJA[::5,::5],units='width',scale=1/0.03,scale_units='width')
axs[3].set_aspect(1.0/axs[3].get_data_ratio()*ratio)
axs[3].set_title('(d) JJA',fontsize=8)

fig.savefig('map_seasonal_wind_ERAi_CMIP5MMM.pdf')
