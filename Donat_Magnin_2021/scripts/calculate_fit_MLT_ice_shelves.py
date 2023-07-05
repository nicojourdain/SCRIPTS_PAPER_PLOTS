import numpy as np
import xarray as xr
import matplotlib
from matplotlib import pyplot as plt

file_MLT='DATA/ERAi_CMIP5anom_ANf/MARv3.9.3_Amundsen_mlt_2088-2117_monthly.nc'
file_T2m='DATA/ERAi_CMIP5anom_ANf/MARv3.9.3_Amundsen_TTz_2088-2117_monthly.nc'
#file_MLT='DATA/ERAi_CMIP5anom_ANf/ICE.2088.2117-mlt.nc'
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
ncA=xr.open_dataset(file_MLT,decode_cf=False)
MLT=ncA['mlt'].values[:,0,:,:]
MLT_isf=np.zeros(np.shape(MLT))

ncB=xr.open_dataset(file_T2m,decode_cf=False)
T2m=ncB['TTz'].values[:,0,:,:]
T2m_isf=np.zeros(np.shape(T2m))

for kmon in np.arange(np.shape(MLT)[0]):
  MLT_isf[kmon,:,:]=MLT[kmon,:,:]*msk_isf
  T2m_isf[kmon,:,:]=T2m[kmon,:,:]*msk_isf
T2m_list=T2m_isf[~np.isnan(MLT_isf)]
MLT_list=MLT_isf[~np.isnan(MLT_isf)]

# fit
# MLT=A.exp(B*T2m) => ln(MLT)=ln(A)+B*T2m
b=np.where((MLT_list>0) & (T2m_list>-10))
T2mb=T2m_list[b]
MLTb=MLT_list[b]
P=np.polyfit(T2mb,np.log(MLTb),1)
print P

fig, ax = plt.subplots()

TT=np.arange(-30,2,0.1)

plt.scatter(T2m_list,MLT_list,s=1,alpha=0.01,color='darkblue')
plt.plot(TT,np.exp(P[1]+P[0]*TT),'--',linewidth=1.0,color='orange')

fig.savefig('fit_MLT.pdf')
