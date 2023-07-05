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

file_pres='DATA/ACCESS_pres_ANg/MARv3.9.3_Amundsen_smb_1989-2009_ACCESS_monthly.nc'
file_futu='DATA/ACCESS_futu_ANh/MARv3.9.3_Amundsen_smb_2080-2100_ACCESS_monthly.nc'
file_anom='DATA/ACCESS_futu_anom_noSICcor_ANi/MARv3.9.3_Amundsen_smb_1989-2009_ANOM_ACCESS_monthly.nc'
file_mona='DATA/ACCESS_futu_anom_withSICcor_ANm/MARv3.9.3_Amundsen_smb_1989-2009_ANOM_ACCESS_monthly.nc'

#====================
ncA=xr.open_dataset(file_pres,decode_cf=False)
smb=ncA['smb'].values[:,0,:,:]
print np.shape(smb)
smbmon=np.zeros((np.shape(smb)[0]))
for kmon in np.arange(np.shape(smb)[0]):
  smbmon[kmon]=np.nansum(np.nansum(smb[kmon,:,:]*msk,axis=1),axis=0)*1.e8*1.e-12*365/12
climsmbA=np.zeros((14))
tmp=np.zeros((12))
for kmon in np.arange(12):
  tmp[kmon] = np.mean(smbmon[kmon::12])
tmp=np.roll(tmp,-4) # to start in May instead of January
climsmbA[0]=tmp[11]
climsmbA[1:13]=tmp
climsmbA[13]=tmp[0]

#====================
ncB=xr.open_dataset(file_futu,decode_cf=False)
smb=ncB['smb'].values[:,0,:,:]
print np.shape(smb)
smbmon=np.zeros((np.shape(smb)[0]))
for kmon in np.arange(np.shape(smb)[0]):
  smbmon[kmon]=np.nansum(np.nansum(smb[kmon,:,:]*msk,axis=1),axis=0)*1.e8*1.e-12*365/12
climsmbB=np.zeros((14))
tmp=np.zeros((12))
for kmon in np.arange(12):
  tmp[kmon] = np.mean(smbmon[kmon::12])
tmp=np.roll(tmp,-4) # to start in May instead of January
climsmbB[0]=tmp[11]
climsmbB[1:13]=tmp
climsmbB[13]=tmp[0]

#====================
ncC=xr.open_dataset(file_anom,decode_cf=False)
smb=ncC['smb'].values[:,0,:,:]
print np.shape(smb)
smbmon=np.zeros((np.shape(smb)[0]))
for kmon in np.arange(np.shape(smb)[0]):
  smbmon[kmon]=np.nansum(np.nansum(smb[kmon,:,:]*msk,axis=1),axis=0)*1.e8*1.e-12*365/12
climsmbC=np.zeros((14))
tmp=np.zeros((12))
for kmon in np.arange(12):
  tmp[kmon] = np.mean(smbmon[kmon::12])
tmp=np.roll(tmp,-4) # to start in May instead of January
climsmbC[0]=tmp[11]
climsmbC[1:13]=tmp
climsmbC[13]=tmp[0]

#====================
ncD=xr.open_dataset(file_mona,decode_cf=False)
smb=ncD['smb'].values[:,0,:,:]
print np.shape(smb)
smbmon=np.zeros((np.shape(smb)[0]))
for kmon in np.arange(np.shape(smb)[0]):
  smbmon[kmon]=np.nansum(np.nansum(smb[kmon,:,:]*msk,axis=1),axis=0)*1.e8*1.e-12*365/12
climsmbD=np.zeros((14))
tmp=np.zeros((12))
for kmon in np.arange(12):
  tmp[kmon] = np.mean(smbmon[kmon::12])
tmp=np.roll(tmp,-4) # to start in May instead of January
climsmbD[0]=tmp[11]
climsmbD[1:13]=tmp
climsmbD[13]=tmp[0]

#====================
fig, ax = plt.subplots()

ax.plot(np.arange(14),climsmbA,label='ACCESS present',linewidth=0.8,color='cornflowerblue')
ax.plot(np.arange(14),climsmbB,label='ACCESS future',linewidth=0.8,color='firebrick')
ax.plot(np.arange(14),climsmbC,'--',label='ACCESS anom. (no sea ice corr.)',linewidth=0.6,color='darkorange')
ax.plot(np.arange(14),climsmbD,label='ACCESS anom. (sea ice corr.)',linewidth=0.8,color='darkorange')
ax.set_xticks(np.arange(1,13,1))
ax.set_xticklabels(['M','J','J','A','S','O','N','D','J','F','M','A'])
plt.xlim(0.5,12.5)
#ax.legend()
ax.set_ylabel('SMB (Gt/month)',size=12)
ax.tick_params(axis='both', which='major', labelsize=10)
ax.text(1.5,50,'(d)',fontsize=18)

print 'annual:'
print np.mean((climsmbD[1:13]-climsmbB[1:13])), np.mean((climsmbC[1:13]-climsmbB[1:13])) # in Gt/month
print np.mean(climsmbB[1:13]-climsmbA[1:13])

error_with_smb_cor_AN=np.mean((climsmbD[1:13]-climsmbB[1:13]))/np.mean((climsmbB[1:13]-climsmbA[1:13]))
error_no_smb_cor_AN=np.mean((climsmbC[1:13]-climsmbB[1:13]))/np.mean((climsmbB[1:13]-climsmbA[1:13]))

print 'relative biases:'
print error_with_smb_cor_AN, error_no_smb_cor_AN

fig.savefig('smb_seasonal_cycle_method_eval.pdf')
