import numpy as np
import xarray as xr

file_pres='DATA/ERAi_ANx/MARv3.9.3_Amundsen_TTz_1979-2017_monthly.nc'
file_futu='DATA/ERAi_CMIP5anom_ANf/MARv3.9.3_Amundsen_TTz_2088-2117_monthly.nc'

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
TTz=ncA['TTz'].values[108:468,0,:,:]  # only 1988-2017 to get the same interannual variability as future projection
TTzmon_ABB=np.zeros((np.shape(TTz)[0]))
TTzmon_COS=np.zeros((np.shape(TTz)[0]))
TTzmon_PIG=np.zeros((np.shape(TTz)[0]))
TTzmon_THW=np.zeros((np.shape(TTz)[0]))
TTzmon_CRO=np.zeros((np.shape(TTz)[0]))
TTzmon_DOT=np.zeros((np.shape(TTz)[0]))
TTzmon_GET=np.zeros((np.shape(TTz)[0]))
for kmon in np.arange(np.shape(TTz)[0]):
  TTzmon_ABB[kmon]=np.nansum(np.nansum(TTz[kmon,:,:]*msk*msk_ABB*(1-msk_isf_ABB),axis=1),axis=0)/np.nansum(np.nansum(msk*msk_ABB*(1-msk_isf_ABB),axis=1),axis=0)
  TTzmon_COS[kmon]=np.nansum(np.nansum(TTz[kmon,:,:]*msk*msk_COS*(1-msk_isf_COS),axis=1),axis=0)/np.nansum(np.nansum(msk*msk_COS*(1-msk_isf_COS),axis=1),axis=0)
  TTzmon_PIG[kmon]=np.nansum(np.nansum(TTz[kmon,:,:]*msk*msk_PIG*(1-msk_isf_PIG),axis=1),axis=0)/np.nansum(np.nansum(msk*msk_PIG*(1-msk_isf_PIG),axis=1),axis=0)
  TTzmon_THW[kmon]=np.nansum(np.nansum(TTz[kmon,:,:]*msk*msk_THW*(1-msk_isf_THW),axis=1),axis=0)/np.nansum(np.nansum(msk*msk_THW*(1-msk_isf_THW),axis=1),axis=0)
  TTzmon_CRO[kmon]=np.nansum(np.nansum(TTz[kmon,:,:]*msk*msk_CRO*(1-msk_isf_CRO),axis=1),axis=0)/np.nansum(np.nansum(msk*msk_CRO*(1-msk_isf_CRO),axis=1),axis=0)
  TTzmon_DOT[kmon]=np.nansum(np.nansum(TTz[kmon,:,:]*msk*msk_DOT*(1-msk_isf_DOT),axis=1),axis=0)/np.nansum(np.nansum(msk*msk_DOT*(1-msk_isf_DOT),axis=1),axis=0)
  TTzmon_GET[kmon]=np.nansum(np.nansum(TTz[kmon,:,:]*msk*msk_GET*(1-msk_isf_GET),axis=1),axis=0)/np.nansum(np.nansum(msk*msk_GET*(1-msk_isf_GET),axis=1),axis=0)
print 'T2m           &       ',np.round(np.mean(TTzmon_ABB),1),'  &       ',\
                               np.round(np.mean(TTzmon_COS),1),'  &       ',\
                               np.round(np.mean(TTzmon_PIG),1),'  &       ',\
                               np.round(np.mean(TTzmon_THW),1),'  &       ',\
                               np.round(np.mean(TTzmon_CRO),1),'  &       ',\
                               np.round(np.mean(TTzmon_DOT),1),'  &       ',\
                               np.round(np.mean(TTzmon_GET),1),' \\\\'
a_ABB=np.mean(TTzmon_ABB)
a_COS=np.mean(TTzmon_COS)
a_PIG=np.mean(TTzmon_PIG)
a_THW=np.mean(TTzmon_THW)
a_CRO=np.mean(TTzmon_CRO)
a_DOT=np.mean(TTzmon_DOT)
a_GET=np.mean(TTzmon_GET)

#====================
ncB=xr.open_dataset(file_futu,decode_cf=False)
TTz=ncB['TTz'].values[:,0,:,:]
TTzmon_ABB=np.zeros((np.shape(TTz)[0]))
TTzmon_COS=np.zeros((np.shape(TTz)[0]))
TTzmon_PIG=np.zeros((np.shape(TTz)[0]))
TTzmon_THW=np.zeros((np.shape(TTz)[0]))
TTzmon_CRO=np.zeros((np.shape(TTz)[0]))
TTzmon_DOT=np.zeros((np.shape(TTz)[0]))
TTzmon_GET=np.zeros((np.shape(TTz)[0]))
for kmon in np.arange(np.shape(TTz)[0]):
  TTzmon_ABB[kmon]=np.nansum(np.nansum(TTz[kmon,:,:]*msk*msk_ABB*(1-msk_isf_ABB),axis=1),axis=0)/np.nansum(np.nansum(msk*msk_ABB*(1-msk_isf_ABB),axis=1),axis=0)
  TTzmon_COS[kmon]=np.nansum(np.nansum(TTz[kmon,:,:]*msk*msk_COS*(1-msk_isf_COS),axis=1),axis=0)/np.nansum(np.nansum(msk*msk_COS*(1-msk_isf_COS),axis=1),axis=0)
  TTzmon_PIG[kmon]=np.nansum(np.nansum(TTz[kmon,:,:]*msk*msk_PIG*(1-msk_isf_PIG),axis=1),axis=0)/np.nansum(np.nansum(msk*msk_PIG*(1-msk_isf_PIG),axis=1),axis=0)
  TTzmon_THW[kmon]=np.nansum(np.nansum(TTz[kmon,:,:]*msk*msk_THW*(1-msk_isf_THW),axis=1),axis=0)/np.nansum(np.nansum(msk*msk_THW*(1-msk_isf_THW),axis=1),axis=0)
  TTzmon_CRO[kmon]=np.nansum(np.nansum(TTz[kmon,:,:]*msk*msk_CRO*(1-msk_isf_CRO),axis=1),axis=0)/np.nansum(np.nansum(msk*msk_CRO*(1-msk_isf_CRO),axis=1),axis=0)
  TTzmon_DOT[kmon]=np.nansum(np.nansum(TTz[kmon,:,:]*msk*msk_DOT*(1-msk_isf_DOT),axis=1),axis=0)/np.nansum(np.nansum(msk*msk_DOT*(1-msk_isf_DOT),axis=1),axis=0)
  TTzmon_GET[kmon]=np.nansum(np.nansum(TTz[kmon,:,:]*msk*msk_GET*(1-msk_isf_GET),axis=1),axis=0)/np.nansum(np.nansum(msk*msk_GET*(1-msk_isf_GET),axis=1),axis=0)
print '              &  {\\bf ',np.round(np.mean(TTzmon_ABB),1),'} & {\\bf ',\
                                np.round(np.mean(TTzmon_COS),1),'} & {\\bf ',\
                                np.round(np.mean(TTzmon_PIG),1),'} & {\\bf ',\
                                np.round(np.mean(TTzmon_THW),1),'} & {\\bf ',\
                                np.round(np.mean(TTzmon_CRO),1),'} & {\\bf ',\
                                np.round(np.mean(TTzmon_DOT),1),'} & {\\bf ',\
                                np.round(np.mean(TTzmon_GET),1),'} \\\\'
b_ABB=np.mean(TTzmon_ABB)
b_COS=np.mean(TTzmon_COS)
b_PIG=np.mean(TTzmon_PIG)
b_THW=np.mean(TTzmon_THW)
b_CRO=np.mean(TTzmon_CRO)
b_DOT=np.mean(TTzmon_DOT)
b_GET=np.mean(TTzmon_GET)

# Delta_T2m, delta_snf/delta_T2m, delta_snf/snf_p/delta_T2m
print 'ABB: ', b_ABB-a_ABB, (50.5-37.0)/(b_ABB-a_ABB), (50.5-37.0)/37.0/(b_ABB-a_ABB)*100.0
print 'COS: ', b_COS-a_COS, (10.0-7.3)/(b_COS-a_COS), (10.0-7.3)/7.3/(b_COS-a_COS)*100.0
print 'PIG: ', b_PIG-a_PIG, (111.3-82.0)/(b_PIG-a_PIG), (111.3-82.0)/82.0/(b_PIG-a_PIG)*100.0
print 'THW: ', b_THW-a_THW, (127.7-95.6)/(b_THW-a_THW), (127.7-95.6)/95.6/(b_THW-a_THW)*100.0
print 'CRO: ', b_CRO-a_CRO, (28.6-21.0)/(b_CRO-a_CRO), (28.6-21.0)/21.0/(b_CRO-a_CRO)*100.0
print 'DOT: ', b_DOT-a_DOT, (23.2-17.3)/(b_DOT-a_DOT), (23.2-17.3)/17.3/(b_DOT-a_DOT)*100.0
print 'GET: ', b_GET-a_GET, (104.5-81.0)/(b_GET-a_GET), (104.5-81.0)/81.0/(b_GET-a_GET)*100.0
