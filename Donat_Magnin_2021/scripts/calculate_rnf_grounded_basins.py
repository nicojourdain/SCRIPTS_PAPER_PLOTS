import numpy as np
import xarray as xr

file_pres='DATA/ERAi_ANx/MARv3.9.3_Amundsen_rnf_1979-2017_monthly.nc'
file_futu='DATA/ERAi_CMIP5anom_ANf/MARv3.9.3_Amundsen_rnf_2088-2117_monthly.nc'

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

#====================
ncA=xr.open_dataset(file_pres,decode_cf=False)
rnf=ncA['rnf'].values[108:468,:,:]  # only 1988-2017 to get the same interannual variability as future projection
rnfmon_ABB=np.zeros((np.shape(rnf)[0]))
rnfmon_COS=np.zeros((np.shape(rnf)[0]))
rnfmon_PIG=np.zeros((np.shape(rnf)[0]))
rnfmon_THW=np.zeros((np.shape(rnf)[0]))
rnfmon_CRO=np.zeros((np.shape(rnf)[0]))
rnfmon_DOT=np.zeros((np.shape(rnf)[0]))
rnfmon_GET=np.zeros((np.shape(rnf)[0]))
for kmon in np.arange(np.shape(rnf)[0]):
  rnfmon_ABB[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_ABB*(1-msk_isf_ABB),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  rnfmon_COS[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_COS*(1-msk_isf_COS),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  rnfmon_PIG[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_PIG*(1-msk_isf_PIG),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  rnfmon_THW[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_THW*(1-msk_isf_THW),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  rnfmon_CRO[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_CRO*(1-msk_isf_CRO),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  rnfmon_DOT[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_DOT*(1-msk_isf_DOT),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  rnfmon_GET[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_GET*(1-msk_isf_GET),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
print 'Rainfall      &       ',np.round(np.mean(rnfmon_ABB),1),'  &       ',\
                               np.round(np.mean(rnfmon_COS),1),'  &       ',\
                               np.round(np.mean(rnfmon_PIG),1),'  &       ',\
                               np.round(np.mean(rnfmon_THW),1),'  &       ',\
                               np.round(np.mean(rnfmon_CRO),1),'  &       ',\
                               np.round(np.mean(rnfmon_DOT),1),'  &       ',\
                               np.round(np.mean(rnfmon_GET),1),' \\\\'

#====================
ncB=xr.open_dataset(file_futu,decode_cf=False)
rnf=ncB['rnf'].values[:,:,:]
rnfmon_ABB=np.zeros((np.shape(rnf)[0]))
rnfmon_COS=np.zeros((np.shape(rnf)[0]))
rnfmon_PIG=np.zeros((np.shape(rnf)[0]))
rnfmon_THW=np.zeros((np.shape(rnf)[0]))
rnfmon_CRO=np.zeros((np.shape(rnf)[0]))
rnfmon_DOT=np.zeros((np.shape(rnf)[0]))
rnfmon_GET=np.zeros((np.shape(rnf)[0]))
for kmon in np.arange(np.shape(rnf)[0]):
  rnfmon_ABB[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_ABB*(1-msk_isf_ABB),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  rnfmon_COS[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_COS*(1-msk_isf_COS),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  rnfmon_PIG[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_PIG*(1-msk_isf_PIG),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  rnfmon_THW[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_THW*(1-msk_isf_THW),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  rnfmon_CRO[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_CRO*(1-msk_isf_CRO),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  rnfmon_DOT[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_DOT*(1-msk_isf_DOT),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
  rnfmon_GET[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_GET*(1-msk_isf_GET),axis=1),axis=0)*365.25*1.e-12*1.e8 # Gt/yr
print '              &  {\\bf ',np.round(np.mean(rnfmon_ABB),1),'} & {\\bf ',\
                                np.round(np.mean(rnfmon_COS),1),'} & {\\bf ',\
                                np.round(np.mean(rnfmon_PIG),1),'} & {\\bf ',\
                                np.round(np.mean(rnfmon_THW),1),'} & {\\bf ',\
                                np.round(np.mean(rnfmon_CRO),1),'} & {\\bf ',\
                                np.round(np.mean(rnfmon_DOT),1),'} & {\\bf ',\
                                np.round(np.mean(rnfmon_GET),1),'} \\\\'
