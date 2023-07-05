import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

file_pres='DATA/ERAi_ANx/MARv3.9.3_Amundsen_LWU_1979-2017_monthly.nc'
file_futu='DATA/ERAi_CMIP5anom_ANf/MARv3.9.3_Amundsen_LWU_2088-2117_monthly.nc'

file_msk='MAR_grid10km.nc'
ncM=xr.open_dataset(file_msk,decode_cf=False)
msk=ncM['MSK'].values[:,:]
msk[msk<100]=np.nan
msk=msk*0.e0+1.e0 # nan over ocean, 1 over ice-sheet (includes ice shelves and grounded ice, but exclude nunatak)

# drainage basins:
msk_ABB=np.load('ANT_MASK_ABBOT.npy'); msk_ABB[ np.isnan(msk_ABB) ]=0
msk_COS=np.load('ANT_MASK_COSGR.npy'); msk_COS[ np.isnan(msk_COS) ]=0
msk_PIG=np.load('ANT_MASK_PINE.npy'); msk_PIG[ np.isnan(msk_PIG) ]=0
msk_THW=np.load('ANT_MASK_THWAIT.npy'); msk_THW[ np.isnan(msk_THW) ]=0
msk_CRO=np.load('ANT_MASK_CROSSON.npy'); msk_CRO[ np.isnan(msk_CRO) ]=0
msk_DOT=np.load('ANT_MASK_DOTSON.npy'); msk_DOT[ np.isnan(msk_DOT) ]=0
msk_GET=np.load('ANT_MASK_GETZ.npy'); msk_GET[ np.isnan(msk_GET) ]=0
msk_basin=msk_ABB+msk_COS+msk_PIG+msk_THW+msk_CRO+msk_DOT+msk_GET

# ice shelves:
msk_isf_ABB=np.load('msk_isf_ABB.npy'); msk_isf_ABB[ np.isnan(msk_isf_ABB) ]=0
msk_isf_COS=np.load('msk_isf_COS.npy'); msk_isf_COS[ np.isnan(msk_isf_COS) ]=0
msk_isf_PIG=np.load('msk_isf_PIG.npy'); msk_isf_PIG[ np.isnan(msk_isf_PIG) ]=0
msk_isf_THW=np.load('msk_isf_THW.npy'); msk_isf_THW[ np.isnan(msk_isf_THW) ]=0
msk_isf_CRO=np.load('msk_isf_CRO.npy'); msk_isf_CRO[ np.isnan(msk_isf_CRO) ]=0
msk_isf_DOT=np.load('msk_isf_DOT.npy'); msk_isf_DOT[ np.isnan(msk_isf_DOT) ]=0
msk_isf_GET=np.load('msk_isf_GET.npy'); msk_isf_GET[ np.isnan(msk_isf_GET) ]=0
msk_isf=msk_isf_ABB+msk_isf_COS+msk_isf_PIG+msk_isf_THW+msk_isf_CRO+msk_isf_DOT+msk_isf_GET

msk=msk*msk_basin*msk_isf

#====================
ncA=xr.open_dataset(file_pres,decode_cf=False)
lwu=ncA['LWU'].values[108:468,:,:]  # only 1988-2017 to get the same interannual variability as future projection
lwumon=np.zeros((np.shape(lwu)[0]))
for kmon in np.arange(np.shape(lwu)[0]):
  lwumon[kmon]=np.nansum(np.nansum(lwu[kmon,:,:]*msk,axis=1),axis=0)/np.nansum(np.nansum(msk,axis=1),axis=0)
climlwuA=np.zeros((14))
tmp=np.zeros((12))
for kmon in np.arange(12):
  tmp[kmon] = np.mean(lwumon[kmon::12])
tmp=np.roll(tmp,-4) # to start in May instead of January
climlwuA[0]=tmp[11]
climlwuA[1:13]=tmp
climlwuA[13]=tmp[0]

#====================
ncB=xr.open_dataset(file_futu,decode_cf=False)
lwu=ncB['LWU'].values[:,:,:]
lwumon=np.zeros((np.shape(lwu)[0]))
for kmon in np.arange(np.shape(lwu)[0]):
  lwumon[kmon]=np.nansum(np.nansum(lwu[kmon,:,:]*msk,axis=1),axis=0)/np.nansum(np.nansum(msk,axis=1),axis=0)
climlwuB=np.zeros((14))
tmp=np.zeros((12))
for kmon in np.arange(12):
  tmp[kmon] = np.mean(lwumon[kmon::12])
tmp=np.roll(tmp,-4) # to start in May instead of January
climlwuB[0]=tmp[11]
climlwuB[1:13]=tmp
climlwuB[13]=tmp[0]

#====================
fig, ax = plt.subplots()

ax.plot(np.arange(14),climlwuA,label='present (ERAinterim)',linewidth=1.0,color='cornflowerblue')
ax.plot(np.arange(14),climlwuB,label='future (ERAinterim + CMIP5 anom)',linewidth=1.0,color='firebrick')
ax.set_xticks(np.arange(1,13,1))
ax.set_xticklabels(['M','J','J','A','S','O','N','D','J','F','M','A'])
plt.xlim(0.5,12.5)
ax.legend()
ax.set_ylabel('Upward longwave (W.m-2)',size=12)
ax.tick_params(axis='both', which='major', labelsize=10)
#ax.text(1,75,'(a)',fontsize=18)

fig.savefig('LWU_seasonal_cycle_ice_shelves_ERAi_CMIP5.pdf')

diff_LWU=climlwuB-climlwuA
np.save('diff_LWU.npy',diff_LWU)
