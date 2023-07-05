import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

file_msk='MAR_grid10km.nc'
ncM=xr.open_dataset(file_msk,decode_cf=False)
msk=ncM['MSK'].values[:,:]
msk[msk<0.1]=np.nan
msk=msk*0.e0+1.e0 # nan over ocean, 1 over ice-sheet

msk_ABB=np.load('ANT_MASK_ABBOT.npy'); msk_ABB[ np.isnan(msk_ABB) ]=0
msk_COS=np.load('ANT_MASK_COSGR.npy'); msk_COS[ np.isnan(msk_COS) ]=0
msk_PIG=np.load('ANT_MASK_PINE.npy'); msk_PIG[ np.isnan(msk_PIG) ]=0
msk_THW=np.load('ANT_MASK_THWAIT.npy'); msk_THW[ np.isnan(msk_THW) ]=0
msk_CRO=np.load('ANT_MASK_CROSSON.npy'); msk_CRO[ np.isnan(msk_CRO) ]=0
msk_DOT=np.load('ANT_MASK_DOTSON.npy'); msk_DOT[ np.isnan(msk_DOT) ]=0
msk_GET=np.load('ANT_MASK_GETZ.npy'); msk_GET[ np.isnan(msk_GET) ]=0
msk_all_drainage_basins=msk_ABB+msk_COS+msk_PIG+msk_THW+msk_CRO+msk_DOT+msk_GET
msk_all_drainage_basins[msk_all_drainage_basins>1]=1
msk_all_drainage_basins[msk_all_drainage_basins < 0.5]=np.nan
msk=msk*msk_all_drainage_basins

file_pres='DATA/ACCESS_pres_ANg/MARv3.9.3_Amundsen_TTz_1989-2009_ACCESS_monthly.nc'
file_futu='DATA/ACCESS_futu_ANh/MARv3.9.3_Amundsen_TTz_2080-2100_ACCESS_monthly.nc'
file_anom='DATA/ACCESS_futu_anom_noSICcor_ANi/MARv3.9.3_Amundsen_TTz_1989-2009_ANOM_ACCESS_monthly.nc'
file_mona='DATA/ACCESS_futu_anom_withSICcor_ANm/MARv3.9.3_Amundsen_TTz_1989-2009_ANOM_ACCESS_monthly.nc'

#====================
ncA=xr.open_dataset(file_pres,decode_cf=False)
t2m=ncA['TTz'].values[:,0,:,:]
t2mmon=np.zeros((np.shape(t2m)[0]))
for kmon in np.arange(np.shape(t2m)[0]):
  t2mmon[kmon]=np.nansum(np.nansum(t2m[kmon,:,:]*msk,axis=1),axis=0)/np.nansum(np.nansum(msk,axis=1),axis=0)
climt2mA=np.zeros((14))
tmp=np.zeros((12))
for kmon in np.arange(12):
  tmp[kmon] = np.mean(t2mmon[kmon::12])
tmp=np.roll(tmp,-4) # to start in May instead of January
climt2mA[0]=tmp[11]
climt2mA[1:13]=tmp
climt2mA[13]=tmp[0]

#====================
ncB=xr.open_dataset(file_futu,decode_cf=False)
t2m=ncB['TTz'].values[:,0,:,:]
t2mmon=np.zeros((np.shape(t2m)[0]))
for kmon in np.arange(np.shape(t2m)[0]):
  t2mmon[kmon]=np.nansum(np.nansum(t2m[kmon,:,:]*msk,axis=1),axis=0)/np.nansum(np.nansum(msk,axis=1),axis=0)
climt2mB=np.zeros((14))
tmp=np.zeros((12))
for kmon in np.arange(12):
  tmp[kmon] = np.mean(t2mmon[kmon::12])
tmp=np.roll(tmp,-4) # to start in May instead of January
climt2mB[0]=tmp[11]
climt2mB[1:13]=tmp
climt2mB[13]=tmp[0]

#====================
ncC=xr.open_dataset(file_anom,decode_cf=False)
t2m=ncC['TTz'].values[:,0,:,:]
t2mmon=np.zeros((np.shape(t2m)[0]))
for kmon in np.arange(np.shape(t2m)[0]):
  t2mmon[kmon]=np.nansum(np.nansum(t2m[kmon,:,:]*msk,axis=1),axis=0)/np.nansum(np.nansum(msk,axis=1),axis=0)
climt2mC=np.zeros((14))
tmp=np.zeros((12))
for kmon in np.arange(12):
  tmp[kmon] = np.mean(t2mmon[kmon::12])
tmp=np.roll(tmp,-4) # to start in May instead of January
climt2mC[0]=tmp[11]
climt2mC[1:13]=tmp
climt2mC[13]=tmp[0]

#====================
ncD=xr.open_dataset(file_mona,decode_cf=False)
t2m=ncD['TTz'].values[:,0,:,:]
t2mmon=np.zeros((np.shape(t2m)[0]))
for kmon in np.arange(np.shape(t2m)[0]):
  t2mmon[kmon]=np.nansum(np.nansum(t2m[kmon,:,:]*msk,axis=1),axis=0)/np.nansum(np.nansum(msk,axis=1),axis=0)
climt2mD=np.zeros((14))
tmp=np.zeros((12))
for kmon in np.arange(12):
  tmp[kmon] = np.mean(t2mmon[kmon::12])
tmp=np.roll(tmp,-4) # to start in May instead of January
climt2mD[0]=tmp[11]
climt2mD[1:13]=tmp
climt2mD[13]=tmp[0]

#====================
fig, ax = plt.subplots()

ax.plot(np.arange(14),climt2mA,label='ACCESS present',linewidth=0.8,color='cornflowerblue')
ax.plot(np.arange(14),climt2mB,label='ACCESS future',linewidth=0.8,color='firebrick')
ax.plot(np.arange(14),climt2mC,'--',label='ACCESS anom. (no sea ice corr.)',linewidth=0.6,color='darkorange')
ax.plot(np.arange(14),climt2mD,label='ACCESS anom. (sea ice corr.)',linewidth=0.8,color='darkorange')
ax.set_xticks(np.arange(1,13,1))
ax.set_xticklabels(['M','J','J','A','S','O','N','D','J','F','M','A'])
plt.xlim(0.5,12.5)
#ax.legend()
ax.set_ylabel('2m air temperature (degC)',size=12)
ax.tick_params(axis='both', which='major', labelsize=10)
ax.text(1,-10,'(b)',fontsize=18)

print 'JJA:'
print np.mean((climt2mD[2:5]-climt2mB[2:5])), np.mean((climt2mC[2:5]-climt2mB[2:5])) # in degC
print np.mean(climt2mB[2:5]-climt2mA[2:5])
print 'DJF:'
print np.mean((climt2mD[8:11]-climt2mB[8:11])), np.mean((climt2mC[8:11]-climt2mB[8:11])) # in degC
print np.mean(climt2mB[8:11]-climt2mA[8:11])

error_with_T2m_cor_JJA=np.mean((climt2mD[2:5]-climt2mB[2:5])/(climt2mB[2:5]-climt2mA[2:5]))
error_no_T2m_cor_JJA=np.mean((climt2mC[2:5]-climt2mB[2:5])/(climt2mB[2:5]-climt2mA[2:5]))

error_with_T2m_cor_DJF=np.mean((climt2mD[8:11]-climt2mB[8:11])/(climt2mB[8:11]-climt2mA[8:11]))
error_no_T2m_cor_DJF=np.mean((climt2mC[8:10]-climt2mB[8:10])/(climt2mB[8:10]-climt2mA[8:10]))

print 'relative biases:'
print error_with_T2m_cor_JJA, error_no_T2m_cor_JJA # in % of climate anomaly
print error_with_T2m_cor_DJF, error_no_T2m_cor_DJF 

fig.savefig('T2m_seasonal_cycle_method_eval.pdf')
