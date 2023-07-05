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

file_pres='DATA/ACCESS_pres_ANg/MARv3.9.3_Amundsen_mlt_1989-2009_ACCESS_monthly.nc'
file_futu='DATA/ACCESS_futu_ANh/MARv3.9.3_Amundsen_mlt_2080-2100_ACCESS_monthly.nc'
file_anom='DATA/ACCESS_futu_anom_noSICcor_ANi/MARv3.9.3_Amundsen_mlt_1989-2009_ANOM_ACCESS_monthly.nc'
file_mona='DATA/ACCESS_futu_anom_withSICcor_ANm/MARv3.9.3_Amundsen_mlt_1989-2009_ANOM_ACCESS_monthly.nc'

#====================
ncA=xr.open_dataset(file_pres,decode_cf=False)
mlt=ncA['mlt'].values[:,0,:,:]
mltmon=np.zeros((np.shape(mlt)[0]))
for kmon in np.arange(np.shape(mlt)[0]):
  mltmon[kmon]=np.nansum(np.nansum(mlt[kmon,:,:]*msk,axis=1),axis=0)*1.e8*1.e-12*365/12
climmltA=np.zeros((14))
tmp=np.zeros((12))
for kmon in np.arange(12):
  tmp[kmon] = np.mean(mltmon[kmon::12])
tmp=np.roll(tmp,-4) # to start in May instead of January
climmltA[0]=tmp[11]
climmltA[1:13]=tmp
climmltA[13]=tmp[0]

#====================
ncB=xr.open_dataset(file_futu,decode_cf=False)
mlt=ncB['mlt'].values[:,0,:,:]
mltmon=np.zeros((np.shape(mlt)[0]))
for kmon in np.arange(np.shape(mlt)[0]):
  mltmon[kmon]=np.nansum(np.nansum(mlt[kmon,:,:]*msk,axis=1),axis=0)*1.e8*1.e-12*365/12
climmltB=np.zeros((14))
tmp=np.zeros((12))
for kmon in np.arange(12):
  tmp[kmon] = np.mean(mltmon[kmon::12])
tmp=np.roll(tmp,-4) # to start in May instead of January
climmltB[0]=tmp[11]
climmltB[1:13]=tmp
climmltB[13]=tmp[0]

#====================
ncC=xr.open_dataset(file_anom,decode_cf=False)
mlt=ncC['mlt'].values[:,0,:,:]
mltmon=np.zeros((np.shape(mlt)[0]))
for kmon in np.arange(np.shape(mlt)[0]):
  mltmon[kmon]=np.nansum(np.nansum(mlt[kmon,:,:]*msk,axis=1),axis=0)*1.e8*1.e-12*365/12
climmltC=np.zeros((14))
tmp=np.zeros((12))
for kmon in np.arange(12):
  tmp[kmon] = np.mean(mltmon[kmon::12])
tmp=np.roll(tmp,-4) # to start in May instead of January
climmltC[0]=tmp[11]
climmltC[1:13]=tmp
climmltC[13]=tmp[0]

#====================
ncD=xr.open_dataset(file_mona,decode_cf=False)
mlt=ncD['mlt'].values[:,0,:,:]
mltmon=np.zeros((np.shape(mlt)[0]))
for kmon in np.arange(np.shape(mlt)[0]):
  mltmon[kmon]=np.nansum(np.nansum(mlt[kmon,:,:]*msk,axis=1),axis=0)*1.e8*1.e-12*365/12
climmltD=np.zeros((14))
tmp=np.zeros((12))
for kmon in np.arange(12):
  tmp[kmon] = np.mean(mltmon[kmon::12])
tmp=np.roll(tmp,-4) # to start in May instead of January
climmltD[0]=tmp[11]
climmltD[1:13]=tmp
climmltD[13]=tmp[0]

#====================
fig, ax = plt.subplots()

ax.plot(np.arange(14),climmltA,label='ACCESS present',linewidth=0.8,color='cornflowerblue')
ax.plot(np.arange(14),climmltB,label='ACCESS future',linewidth=0.8,color='firebrick')
ax.plot(np.arange(14),climmltC,'--',label='ACCESS anom. (no sea ice corr.)',linewidth=0.6,color='darkorange')
ax.plot(np.arange(14),climmltD,label='ACCESS anom. (sea ice corr.)',linewidth=0.8,color='darkorange')
ax.set_xticks(np.arange(1,13,1))
ax.set_xticklabels(['M','J','J','A','S','O','N','D','J','F','M','A'])
plt.xlim(0.5,12.5)
#ax.legend()
ax.set_ylabel('Surface melting (Gt/month)',size=12)
ax.tick_params(axis='both', which='major', labelsize=10)
ax.text(1,25,'(c)',fontsize=18)

print 'DJF:'
print np.mean((climmltD[8:11]-climmltB[8:11])), np.mean((climmltC[8:11]-climmltB[8:11])) # in degC
print np.mean(climmltB[8:11]-climmltA[8:11])
print 'January (peak):'
print np.mean((climmltD[9]-climmltB[9])), np.mean((climmltC[9]-climmltB[9])) # in degC
print np.mean(climmltB[9]-climmltA[9])

error_with_mlt_cor_DJF=np.mean((climmltD[8:11]-climmltB[8:11])/(climmltB[8:11]-climmltA[8:11]))
error_no_mlt_cor_DJF=np.mean((climmltC[8:11]-climmltB[8:11])/(climmltB[8:11]-climmltA[8:11]))
error_with_mlt_cor_Jan=np.mean((climmltD[9]-climmltB[9])/(climmltB[9]-climmltA[9]))
error_no_mlt_cor_Jan=np.mean((climmltC[9]-climmltB[9])/(climmltB[9]-climmltA[9]))

print 'relative biases:'
print error_with_mlt_cor_DJF, error_no_mlt_cor_DJF
print error_with_mlt_cor_Jan, error_no_mlt_cor_Jan

fig.savefig('mlt_seasonal_cycle_method_eval.pdf')
