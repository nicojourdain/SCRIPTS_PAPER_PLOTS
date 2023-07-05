import numpy as np
import xarray as xr

file_pres='DATA/ERAi_ANx/MARv3.9.3_Amundsen_smb_1979-2017_monthly.nc'
file_futu='DATA/ERAi_CMIP5anom_ANf/MARv3.9.3_Amundsen_smb_2088-2117_monthly.nc'

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

# ice shelves:
msk_isf_ABB=np.load('msk_isf_ABB.npy'); msk_isf_ABB[ np.isnan(msk_isf_ABB) ]=0
msk_isf_COS=np.load('msk_isf_COS.npy'); msk_isf_COS[ np.isnan(msk_isf_COS) ]=0
msk_isf_PIG=np.load('msk_isf_PIG.npy'); msk_isf_PIG[ np.isnan(msk_isf_PIG) ]=0
msk_isf_THW=np.load('msk_isf_THW.npy'); msk_isf_THW[ np.isnan(msk_isf_THW) ]=0
msk_isf_CRO=np.load('msk_isf_CRO.npy'); msk_isf_CRO[ np.isnan(msk_isf_CRO) ]=0
msk_isf_DOT=np.load('msk_isf_DOT.npy'); msk_isf_DOT[ np.isnan(msk_isf_DOT) ]=0
msk_isf_GET=np.load('msk_isf_GET.npy'); msk_isf_GET[ np.isnan(msk_isf_GET) ]=0


print('areas (10^3 km^2):')
print('ABB:',np.nansum(np.nansum(msk*msk_ABB*(1-msk_isf_ABB),axis=1),axis=0)*1.e-1)
print('COS:',np.nansum(np.nansum(msk*msk_COS*(1-msk_isf_COS),axis=1),axis=0)*1.e-1)
print('PIG:',np.nansum(np.nansum(msk*msk_PIG*(1-msk_isf_PIG),axis=1),axis=0)*1.e-1)
print('THW:',np.nansum(np.nansum(msk*msk_THW*(1-msk_isf_THW),axis=1),axis=0)*1.e-1)
print('CRO:',np.nansum(np.nansum(msk*msk_CRO*(1-msk_isf_CRO),axis=1),axis=0)*1.e-1)
print('DOT:',np.nansum(np.nansum(msk*msk_DOT*(1-msk_isf_DOT),axis=1),axis=0)*1.e-1)
print('GET:',np.nansum(np.nansum(msk*msk_GET*(1-msk_isf_GET),axis=1),axis=0)*1.e-1)
print('-----------------------------------')

#====================
ncA=xr.open_dataset(file_pres,decode_cf=False)
smb=ncA['smb'].values[108:468,0,:,:]  # only 1988-2017 to get the same interannual variability as future projection
smbmon_ABB=np.zeros((np.shape(smb)[0]))
smbmon_COS=np.zeros((np.shape(smb)[0]))
smbmon_PIG=np.zeros((np.shape(smb)[0]))
smbmon_THW=np.zeros((np.shape(smb)[0]))
smbmon_CRO=np.zeros((np.shape(smb)[0]))
smbmon_DOT=np.zeros((np.shape(smb)[0]))
smbmon_GET=np.zeros((np.shape(smb)[0]))
for kmon in np.arange(np.shape(smb)[0]):
  smbmon_ABB[kmon]=np.nansum(np.nansum(smb[kmon,:,:]*msk*msk_ABB*(1-msk_isf_ABB),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  smbmon_COS[kmon]=np.nansum(np.nansum(smb[kmon,:,:]*msk*msk_COS*(1-msk_isf_COS),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  smbmon_PIG[kmon]=np.nansum(np.nansum(smb[kmon,:,:]*msk*msk_PIG*(1-msk_isf_PIG),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  smbmon_THW[kmon]=np.nansum(np.nansum(smb[kmon,:,:]*msk*msk_THW*(1-msk_isf_THW),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  smbmon_CRO[kmon]=np.nansum(np.nansum(smb[kmon,:,:]*msk*msk_CRO*(1-msk_isf_CRO),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  smbmon_DOT[kmon]=np.nansum(np.nansum(smb[kmon,:,:]*msk*msk_DOT*(1-msk_isf_DOT),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  smbmon_GET[kmon]=np.nansum(np.nansum(smb[kmon,:,:]*msk*msk_GET*(1-msk_isf_GET),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
print(np.mean(smbmon_ABB)+np.mean(smbmon_COS)+np.mean(smbmon_PIG)+np.mean(smbmon_THW)+np.mean(smbmon_CRO)+np.mean(smbmon_DOT)+np.mean(smbmon_GET))
print('SMB           &       ',np.round(np.mean(smbmon_ABB),1),'  &       ',\
                               np.round(np.mean(smbmon_COS),1),'  &       ',\
                               np.round(np.mean(smbmon_PIG),1),'  &       ',\
                               np.round(np.mean(smbmon_THW),1),'  &       ',\
                               np.round(np.mean(smbmon_CRO),1),'  &       ',\
                               np.round(np.mean(smbmon_DOT),1),'  &       ',\
                               np.round(np.mean(smbmon_GET),1),' \\\\')

#====================
ncB=xr.open_dataset(file_futu,decode_cf=False)
smb=ncB['smb'].values[:,0,:,:]
smbmon_ABB=np.zeros((np.shape(smb)[0]))
smbmon_COS=np.zeros((np.shape(smb)[0]))
smbmon_PIG=np.zeros((np.shape(smb)[0]))
smbmon_THW=np.zeros((np.shape(smb)[0]))
smbmon_CRO=np.zeros((np.shape(smb)[0]))
smbmon_DOT=np.zeros((np.shape(smb)[0]))
smbmon_GET=np.zeros((np.shape(smb)[0]))
for kmon in np.arange(np.shape(smb)[0]):
  smbmon_ABB[kmon]=np.nansum(np.nansum(smb[kmon,:,:]*msk*msk_ABB*(1-msk_isf_ABB),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  smbmon_COS[kmon]=np.nansum(np.nansum(smb[kmon,:,:]*msk*msk_COS*(1-msk_isf_COS),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  smbmon_PIG[kmon]=np.nansum(np.nansum(smb[kmon,:,:]*msk*msk_PIG*(1-msk_isf_PIG),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  smbmon_THW[kmon]=np.nansum(np.nansum(smb[kmon,:,:]*msk*msk_THW*(1-msk_isf_THW),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  smbmon_CRO[kmon]=np.nansum(np.nansum(smb[kmon,:,:]*msk*msk_CRO*(1-msk_isf_CRO),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  smbmon_DOT[kmon]=np.nansum(np.nansum(smb[kmon,:,:]*msk*msk_DOT*(1-msk_isf_DOT),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  smbmon_GET[kmon]=np.nansum(np.nansum(smb[kmon,:,:]*msk*msk_GET*(1-msk_isf_GET),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
print('              &  {\\bf ',np.round(np.mean(smbmon_ABB),1),'} & {\\bf ',\
                                np.round(np.mean(smbmon_COS),1),'} & {\\bf ',\
                                np.round(np.mean(smbmon_PIG),1),'} & {\\bf ',\
                                np.round(np.mean(smbmon_THW),1),'} & {\\bf ',\
                                np.round(np.mean(smbmon_CRO),1),'} & {\\bf ',\
                                np.round(np.mean(smbmon_DOT),1),'} & {\\bf ',\
                                np.round(np.mean(smbmon_GET),1),'} \\\\')
print(np.mean(smbmon_ABB)+np.mean(smbmon_COS)+np.mean(smbmon_PIG)+np.mean(smbmon_THW)+np.mean(smbmon_CRO)+np.mean(smbmon_DOT)+np.mean(smbmon_GET))
