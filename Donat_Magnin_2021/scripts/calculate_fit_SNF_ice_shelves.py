import numpy as np
import xarray as xr
import matplotlib
from matplotlib import pyplot as plt

file_SNF='DATA/ERAi_CMIP5anom_ANf/MARv3.9.3_Amundsen_snf_2088-2117_monthly.nc'
file_T2m='DATA/ERAi_CMIP5anom_ANf/MARv3.9.3_Amundsen_TTz_2088-2117_monthly.nc'
#file_SNF='DATA/ERAi_CMIP5anom_ANf/ICE.2088.2117-snf.nc'
#file_T2m='DATA/ERAi_CMIP5anom_ANf/ICE.2088.2117-TTz.nc'


# ice shelves:
msk_isf_ABB=np.load('msk_isf_ABB.npy'); msk_isf_ABB[ np.isnan(msk_isf_ABB) ]=0
msk_isf_COS=np.load('msk_isf_COS.npy'); msk_isf_COS[ np.isnan(msk_isf_COS) ]=0
msk_isf_PIG=np.load('msk_isf_PIG.npy'); msk_isf_PIG[ np.isnan(msk_isf_PIG) ]=0
msk_isf_THW=np.load('msk_isf_THW.npy'); msk_isf_THW[ np.isnan(msk_isf_THW) ]=0
msk_isf_CRO=np.load('msk_isf_CRO.npy'); msk_isf_CRO[ np.isnan(msk_isf_CRO) ]=0
msk_isf_DOT=np.load('msk_isf_DOT.npy'); msk_isf_DOT[ np.isnan(msk_isf_DOT) ]=0
msk_isf_GET=np.load('msk_isf_GET.npy'); msk_isf_GET[ np.isnan(msk_isf_GET) ]=0
msk_isf=msk_isf_ABB+msk_isf_COS+msk_isf_PIG+msk_isf_THW+msk_isf_CRO+msk_isf_DOT+msk_isf_GET
msk_isf[msk_isf==0]=np.nan

#====================
ncA=xr.open_dataset(file_SNF,decode_cf=False)
SNF=ncA['snf'].values[:,:,:]
SNF_isf=np.zeros(np.shape(SNF))

ncB=xr.open_dataset(file_T2m,decode_cf=False)
T2m=ncB['TTz'].values[:,0,:,:]
T2m_isf=np.zeros(np.shape(T2m))

for kmon in np.arange(np.shape(SNF)[0]):
  SNF_isf[kmon,:,:]=SNF[kmon,:,:]*msk_isf
  T2m_isf[kmon,:,:]=T2m[kmon,:,:]*msk_isf
T2m_list=T2m_isf[~np.isnan(SNF_isf)]
SNF_list=SNF_isf[~np.isnan(SNF_isf)]

fig, ax = plt.subplots()

plt.scatter(T2m_list,SNF_list,s=2)

fig.savefig('fit_SNF.pdf')
