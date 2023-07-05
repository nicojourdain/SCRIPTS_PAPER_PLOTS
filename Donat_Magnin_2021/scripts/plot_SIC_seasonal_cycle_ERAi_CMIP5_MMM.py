import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

file_msk='MAR_grid10km.nc'

file_pres='DATA/ERAi_ANx/MARv3.9.3_Amundsen_SIC_1979-2017_monthly.nc'
file_futu='DATA/ERAi_CMIP5anom_ANf/MARv3.9.3_Amundsen_SIC_2088-2117_monthly.nc'

ncM=xr.open_dataset(file_msk,decode_cf=False)
lat=ncM['LAT'].values[:,:]
msk=ncM['MSK'].values[:,:]
msk[msk>0.01]=np.nan
msk=msk*0.e0+1.e0
msklat=msk*1.e0
msklat[lat>-70.0]=np.nan

#====================
ncA=xr.open_dataset(file_pres,decode_cf=False)
sic=ncA['SIC'].values[108:468,0,:,:]  # only 1988-2017 to get the same interannual variability as future projection
sicmon=np.zeros((np.shape(sic)[0]))
sicmon70=np.zeros((np.shape(sic)[0]))
for kmon in np.arange(np.shape(sic)[0]):
  sicmon[kmon]=np.nansum(np.nansum(sic[kmon,:,:]*msk,axis=1),axis=0)/np.nansum(np.nansum(msk,axis=1),axis=0)
  sicmon70[kmon]=np.nansum(np.nansum(sic[kmon,:,:]*msklat,axis=1),axis=0)/np.nansum(np.nansum(msklat,axis=1),axis=0)
climsicA=np.zeros((14))
climsicA70=np.zeros((14))
tmp=np.zeros((12))
tmp70=np.zeros((12))
for kmon in np.arange(12):
  tmp[kmon] = np.mean(sicmon[kmon::12])
  tmp70[kmon] = np.mean(sicmon70[kmon::12])
tmp=np.roll(tmp,-4) # to start in May instead of January
tmp70=np.roll(tmp70,-4) # to start in May instead of January
climsicA[0]=tmp[11]
climsicA[1:13]=tmp
climsicA[13]=tmp[0]
climsicA70[0]=tmp70[11]
climsicA70[1:13]=tmp70
climsicA70[13]=tmp70[0]

#====================
ncB=xr.open_dataset(file_futu,decode_cf=False)
sic=ncB['SIC'].values[:,0,:,:]
sicmon=np.zeros((np.shape(sic)[0]))
sicmon70=np.zeros((np.shape(sic)[0]))
for kmon in np.arange(np.shape(sic)[0]):
  sicmon[kmon]=np.nansum(np.nansum(sic[kmon,:,:]*msk,axis=1),axis=0)/np.nansum(np.nansum(msk,axis=1),axis=0)
  sicmon70[kmon]=np.nansum(np.nansum(sic[kmon,:,:]*msklat,axis=1),axis=0)/np.nansum(np.nansum(msklat,axis=1),axis=0)
climsicB=np.zeros((14))
climsicB70=np.zeros((14))
tmp=np.zeros((12))
tmp70=np.zeros((12))
for kmon in np.arange(12):
  tmp[kmon] = np.mean(sicmon[kmon::12])
  tmp70[kmon] = np.mean(sicmon70[kmon::12])
tmp=np.roll(tmp,-4) # to start in May instead of January
tmp70=np.roll(tmp70,-4) # to start in May instead of January
climsicB[0]=tmp[11]
climsicB[1:13]=tmp
climsicB[13]=tmp[0]
climsicB70[0]=tmp70[11]
climsicB70[1:13]=tmp70
climsicB70[13]=tmp70[0]

#====================
fig, ax = plt.subplots()

ax.plot(np.arange(14),climsicA,label='present (ERAinterim)',linewidth=1.0,color='cornflowerblue')
ax.plot(np.arange(14),climsicB,label='future (ERAinterim + CMIP5 anom)',linewidth=1.0,color='firebrick')
ax.plot(np.arange(14),climsicA70,'--',label='present, southward of 70S',linewidth=0.75,color='cornflowerblue')
ax.plot(np.arange(14),climsicB70,'--',label='future, southward of 70S',linewidth=0.75,color='firebrick')
ax.set_xticks(np.arange(1,13,1))
ax.set_xticklabels(['M','J','J','A','S','O','N','D','J','F','M','A'])
plt.xlim(0.5,12.5)
ax.legend()
ax.set_ylabel('sea ice cover (%)',size=12)
ax.tick_params(axis='both', which='major', labelsize=10)
#ax.text(1,75,'(a)',fontsize=18)

fig.savefig('SIC_seasonal_cycle_ERAi_CMIP5.pdf')
