import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

file_msk='MAR_grid10km.nc'

file_pres='DATA/ACCESS_pres_ANg/MARv3.9.3_Amundsen_SIC_1989-2009_ACCESS_monthly.nc'
file_futu='DATA/ACCESS_futu_ANh/MARv3.9.3_Amundsen_SIC_2080-2100_ACCESS_monthly.nc'
file_anom='DATA/ACCESS_futu_anom_noSICcor_ANi/MARv3.9.3_Amundsen_SIC_1989-2009_ANOM_ACCESS_monthly.nc'
file_mona='DATA/ACCESS_futu_anom_withSICcor_ANm/MARv3.9.3_Amundsen_SIC_1989-2009_ANOM_ACCESS_monthly.nc'

ncM=xr.open_dataset(file_msk,decode_cf=False)
msk=ncM['MSK'].values[:,:]
msk[msk>0.01]=np.nan
msk=msk*0.e0+1.e0

#====================
ncA=xr.open_dataset(file_pres,decode_cf=False)
sic=ncA['SIC'].values[:,0,:,:]
sicmon=np.zeros((np.shape(sic)[0]))
for kmon in np.arange(np.shape(sic)[0]):
  sicmon[kmon]=np.nansum(np.nansum(sic[kmon,:,:]*msk,axis=1),axis=0)/np.nansum(np.nansum(msk,axis=1),axis=0)
climsicA=np.zeros((14))
tmp=np.zeros((12))
for kmon in np.arange(12):
  tmp[kmon] = np.mean(sicmon[kmon::12])
tmp=np.roll(tmp,-4) # to start in May instead of January
climsicA[0]=tmp[11]
climsicA[1:13]=tmp
climsicA[13]=tmp[0]

#====================
ncB=xr.open_dataset(file_futu,decode_cf=False)
sic=ncB['SIC'].values[:,0,:,:]
sicmon=np.zeros((np.shape(sic)[0]))
for kmon in np.arange(np.shape(sic)[0]):
  sicmon[kmon]=np.nansum(np.nansum(sic[kmon,:,:]*msk,axis=1),axis=0)/np.nansum(np.nansum(msk,axis=1),axis=0)
climsicB=np.zeros((14))
tmp=np.zeros((12))
for kmon in np.arange(12):
  tmp[kmon] = np.mean(sicmon[kmon::12])
tmp=np.roll(tmp,-4) # to start in May instead of January
climsicB[0]=tmp[11]
climsicB[1:13]=tmp
climsicB[13]=tmp[0]

#====================
ncC=xr.open_dataset(file_anom,decode_cf=False)
sic=ncC['SIC'].values[:,0,:,:]
sicmon=np.zeros((np.shape(sic)[0]))
for kmon in np.arange(np.shape(sic)[0]):
  sicmon[kmon]=np.nansum(np.nansum(sic[kmon,:,:]*msk,axis=1),axis=0)/np.nansum(np.nansum(msk,axis=1),axis=0)
climsicC=np.zeros((14))
tmp=np.zeros((12))
for kmon in np.arange(12):
  tmp[kmon] = np.mean(sicmon[kmon::12])
tmp=np.roll(tmp,-4) # to start in May instead of January
climsicC[0]=tmp[11]
climsicC[1:13]=tmp
climsicC[13]=tmp[0]

#====================
ncD=xr.open_dataset(file_mona,decode_cf=False)
sic=ncD['SIC'].values[:,0,:,:]
sicmon=np.zeros((np.shape(sic)[0]))
for kmon in np.arange(np.shape(sic)[0]):
  sicmon[kmon]=np.nansum(np.nansum(sic[kmon,:,:]*msk,axis=1),axis=0)/np.nansum(np.nansum(msk,axis=1),axis=0)
climsicD=np.zeros((14))
tmp=np.zeros((12))
for kmon in np.arange(12):
  tmp[kmon] = np.mean(sicmon[kmon::12])
tmp=np.roll(tmp,-4) # to start in May instead of January
climsicD[0]=tmp[11]
climsicD[1:13]=tmp
climsicD[13]=tmp[0]

#====================
fig, ax = plt.subplots()

ax.plot(np.arange(14),climsicA,label='ACCESS present',linewidth=0.8,color='cornflowerblue')
ax.plot(np.arange(14),climsicB,label='ACCESS future',linewidth=0.8,color='firebrick')
ax.plot(np.arange(14),climsicC,'--',label='ACCESS anom. (no sea ice corr.)',linewidth=0.6,color='darkorange')
ax.plot(np.arange(14),climsicD,label='ACCESS anom. (sea ice corr.)',linewidth=0.8,color='darkorange')
ax.set_xticks(np.arange(1,13,1))
ax.set_xticklabels(['M','J','J','A','S','O','N','D','J','F','M','A'])
plt.xlim(0.5,12.5)
ax.legend()
ax.set_ylabel('sea ice cover (%)',size=12)
ax.tick_params(axis='both', which='major', labelsize=10)
ax.text(1,75,'(a)',fontsize=18)

error_with_SIC_cor_SON=np.mean((climsicD[5:8]-climsicB[5:8])/(climsicB[5:8]-climsicA[5:8]))
error_no_SIC_cor_SON=np.mean((climsicC[5:8]-climsicB[5:8])/(climsicB[5:8]-climsicA[5:8]))

error_with_SIC_cor_DJF=np.mean((climsicD[8:11]-climsicB[8:11])/(climsicB[8:11]-climsicA[8:11]))
error_no_SIC_cor_DJF=np.mean((climsicC[8:11]-climsicB[8:11])/(climsicB[8:11]-climsicA[8:11]))

print error_with_SIC_cor_SON, error_no_SIC_cor_SON
print error_with_SIC_cor_DJF, error_no_SIC_cor_DJF

fig.savefig('SIC_seasonal_cycle_method_eval.pdf')
