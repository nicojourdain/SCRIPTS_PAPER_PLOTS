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
  rnfmon_ABB[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_ABB*msk_isf_ABB,axis=1),axis=0)/np.nansum(np.nansum(msk*msk_ABB*msk_isf_ABB,axis=1),axis=0)*365.25 # mm.w.eq/yr
  rnfmon_COS[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_COS*msk_isf_COS,axis=1),axis=0)/np.nansum(np.nansum(msk*msk_COS*msk_isf_COS,axis=1),axis=0)*365.25
  rnfmon_PIG[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_PIG*msk_isf_PIG,axis=1),axis=0)/np.nansum(np.nansum(msk*msk_PIG*msk_isf_PIG,axis=1),axis=0)*365.25
  rnfmon_THW[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_THW*msk_isf_THW,axis=1),axis=0)/np.nansum(np.nansum(msk*msk_THW*msk_isf_THW,axis=1),axis=0)*365.25
  rnfmon_CRO[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_CRO*msk_isf_CRO,axis=1),axis=0)/np.nansum(np.nansum(msk*msk_CRO*msk_isf_CRO,axis=1),axis=0)*365.25
  rnfmon_DOT[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_DOT*msk_isf_DOT,axis=1),axis=0)/np.nansum(np.nansum(msk*msk_DOT*msk_isf_DOT,axis=1),axis=0)*365.25
  rnfmon_GET[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_GET*msk_isf_GET,axis=1),axis=0)/np.nansum(np.nansum(msk*msk_GET*msk_isf_GET,axis=1),axis=0)*365.25
print 'Rainfall      &       ',np.round(np.mean(rnfmon_ABB),0),'  &       ',\
                               np.round(np.mean(rnfmon_COS),0),'  &       ',\
                               np.round(np.mean(rnfmon_PIG),0),'  &       ',\
                               np.round(np.mean(rnfmon_THW),0),'  &       ',\
                               np.round(np.mean(rnfmon_CRO),0),'  &       ',\
                               np.round(np.mean(rnfmon_DOT),0),'  &       ',\
                               np.round(np.mean(rnfmon_GET),0),' \\\\'

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
  rnfmon_ABB[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_ABB*msk_isf_ABB,axis=1),axis=0)/np.nansum(np.nansum(msk*msk_ABB*msk_isf_ABB,axis=1),axis=0)*365.25 # mm.w.eq/yr
  rnfmon_COS[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_COS*msk_isf_COS,axis=1),axis=0)/np.nansum(np.nansum(msk*msk_COS*msk_isf_COS,axis=1),axis=0)*365.25
  rnfmon_PIG[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_PIG*msk_isf_PIG,axis=1),axis=0)/np.nansum(np.nansum(msk*msk_PIG*msk_isf_PIG,axis=1),axis=0)*365.25
  rnfmon_THW[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_THW*msk_isf_THW,axis=1),axis=0)/np.nansum(np.nansum(msk*msk_THW*msk_isf_THW,axis=1),axis=0)*365.25
  rnfmon_CRO[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_CRO*msk_isf_CRO,axis=1),axis=0)/np.nansum(np.nansum(msk*msk_CRO*msk_isf_CRO,axis=1),axis=0)*365.25
  rnfmon_DOT[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_DOT*msk_isf_DOT,axis=1),axis=0)/np.nansum(np.nansum(msk*msk_DOT*msk_isf_DOT,axis=1),axis=0)*365.25
  rnfmon_GET[kmon]=np.nansum(np.nansum(rnf[kmon,:,:]*msk*msk_GET*msk_isf_GET,axis=1),axis=0)/np.nansum(np.nansum(msk*msk_GET*msk_isf_GET,axis=1),axis=0)*365.25
print '              &  {\\bf ',np.round(np.mean(rnfmon_ABB),0),'} & {\\bf ',\
                                np.round(np.mean(rnfmon_COS),0),'} & {\\bf ',\
                                np.round(np.mean(rnfmon_PIG),0),'} & {\\bf ',\
                                np.round(np.mean(rnfmon_THW),0),'} & {\\bf ',\
                                np.round(np.mean(rnfmon_CRO),0),'} & {\\bf ',\
                                np.round(np.mean(rnfmon_DOT),0),'} & {\\bf ',\
                                np.round(np.mean(rnfmon_GET),0),'} \\\\'
