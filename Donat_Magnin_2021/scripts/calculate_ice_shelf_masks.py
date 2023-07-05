import numpy as np
import xarray as xr

msk_ABB=np.load('ANT_MASK_ABBOT.npy')
msk_COS=np.load('ANT_MASK_COSGR.npy')
msk_PIG=np.load('ANT_MASK_PINE.npy')
msk_THW=np.load('ANT_MASK_THWAIT.npy')
msk_CRO=np.load('ANT_MASK_CROSSON.npy')
msk_DOT=np.load('ANT_MASK_DOTSON.npy')
msk_GET=np.load('ANT_MASK_GETZ.npy')

ncG=xr.open_dataset('MAR_grid10km.nc')
sh=ncG['SH'].values[:,:]

msk_isf=sh*1
msk_isf[sh<10.0]=np.nan
msk_isf[sh>100.0]=np.nan
msk_isf[msk_isf>1]=1

print np.shape(msk_GET), np.nanmax(msk_GET), np.nanmin(msk_GET)
print np.shape(msk_isf), np.nanmax(msk_isf), np.nanmin(msk_isf)

msk_isf_ABB=msk_ABB*msk_isf
msk_isf_COS=msk_COS*msk_isf
msk_isf_PIG=msk_PIG*msk_isf
msk_isf_THW=msk_THW*msk_isf
msk_isf_CRO=msk_CRO*msk_isf
msk_isf_DOT=msk_DOT*msk_isf
msk_isf_GET=msk_GET*msk_isf

np.save('msk_isf_ABB.npy',msk_isf_ABB)
np.save('msk_isf_COS.npy',msk_isf_COS)
np.save('msk_isf_PIG.npy',msk_isf_PIG)
np.save('msk_isf_THW.npy',msk_isf_THW)
np.save('msk_isf_CRO.npy',msk_isf_CRO)
np.save('msk_isf_DOT.npy',msk_isf_DOT)
np.save('msk_isf_GET.npy',msk_isf_GET)
