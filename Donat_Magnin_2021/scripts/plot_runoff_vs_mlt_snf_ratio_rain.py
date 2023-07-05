import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

file_msk='MAR_grid10km.nc'
ncM=xr.open_dataset(file_msk,decode_cf=False)
msk=ncM['MSK'].values[:,:]
msk[msk<0.9]=np.nan
msk=msk*0.e0+1.e0 # nan over ocean, 1 over ice-sheet

#msk_ABB=np.load('ANT_MASK_ABBOT.npy'); msk_ABB[ np.isnan(msk_ABB) ]=0
msk_COS=np.load('ANT_MASK_COSGR.npy'); msk_COS[ np.isnan(msk_COS) ]=0
msk_PIG=np.load('ANT_MASK_PINE.npy'); msk_PIG[ np.isnan(msk_PIG) ]=0
msk_THW=np.load('ANT_MASK_THWAIT.npy'); msk_THW[ np.isnan(msk_THW) ]=0
msk_CRO=np.load('ANT_MASK_CROSSON.npy'); msk_CRO[ np.isnan(msk_CRO) ]=0
msk_DOT=np.load('ANT_MASK_DOTSON.npy'); msk_DOT[ np.isnan(msk_DOT) ]=0
msk_GET=np.load('ANT_MASK_GETZ.npy'); msk_GET[ np.isnan(msk_GET) ]=0
#msk_all_drainage_basins=msk_ABB+msk_COS+msk_PIG+msk_THW+msk_CRO+msk_DOT+msk_GET
msk_all_drainage_basins=msk_COS+msk_PIG+msk_THW+msk_CRO+msk_DOT+msk_GET  ####### ABBOT NOT CONSIDERED HERE #######
msk_all_drainage_basins[msk_all_drainage_basins>1]=1
msk_all_drainage_basins[msk_all_drainage_basins < 0.5]=np.nan
msk=msk*msk_all_drainage_basins

pref_pres='DATA/ERAi_ANx/MARv3.9.3_Amundsen_'
suff_pres='_1979-2017_monthly.nc'
pref_futu='DATA/ERAi_CMIP5anom_ANf/MARv3.9.3_Amundsen_'
suff_futu='_2088-2117_monthly.nc'

##====================
#file_mlt=pref_pres+'mlt'+suff_pres
#ncA1=xr.open_dataset(file_mlt,decode_cf=False)
#mlt=ncA1['mlt'].values[:,0,:,:]
#mlt_clim=np.reshape( np.nanmean(mlt,axis=0)*msk, np.shape(mlt)[1]*np.shape(mlt)[2] ) 
#
#file_rof=pref_pres+'rof'+suff_pres
#ncA2=xr.open_dataset(file_rof,decode_cf=False)
#rof=ncA2['rof'].values[:,0,:,:]
#rof_clim=np.reshape( np.nanmean(rof,axis=0)*msk*365*1.e-3, np.shape(rof)[1]*np.shape(rof)[2] ) # in m.w.e/yr
#
#file_snf=pref_pres+'snf'+suff_pres
#ncA3=xr.open_dataset(file_snf,decode_cf=False)
#snf=ncA3['snf'].values[:,:,:]
#snf_clim=np.reshape( np.nanmean(snf,axis=0)*msk, np.shape(snf)[1]*np.shape(snf)[2] )
#
#tmp=mlt_clim*rof_clim*snf_clim
#mltA=mlt_clim[ ~np.isnan(tmp) ]
#rofA=rof_clim[ ~np.isnan(tmp) ]
#snfA=snf_clim[ ~np.isnan(tmp) ]

#====================
file_mlt=pref_futu+'mlt'+suff_futu
ncB1=xr.open_dataset(file_mlt,decode_cf=False)
mlt=ncB1['mlt'].values[:,0,:,:]
mlt_clim=np.reshape( np.nanmean(mlt,axis=0)*msk, np.shape(mlt)[1]*np.shape(mlt)[2] ) 

file_rof=pref_futu+'rof'+suff_futu
ncB2=xr.open_dataset(file_rof,decode_cf=False)
rof=ncB2['rof'].values[:,0,:,:]
rof_clim=np.reshape( np.nanmean(rof,axis=0)*msk*365*1.e-3, np.shape(rof)[1]*np.shape(rof)[2] ) # in m.w.e/yr

file_snf=pref_futu+'snf'+suff_futu
ncB3=xr.open_dataset(file_snf,decode_cf=False)
snf=ncB3['snf'].values[:,:,:]
snf_clim=np.reshape( np.nanmean(snf,axis=0)*msk, np.shape(snf)[1]*np.shape(snf)[2] )

file_rnf=pref_futu+'rnf'+suff_futu
ncB4=xr.open_dataset(file_rnf,decode_cf=False)
rnf=ncB4['rnf'].values[:,:,:]
rnf_clim=np.reshape( np.nanmean(rnf,axis=0)*msk, np.shape(rnf)[1]*np.shape(rnf)[2] )

tmp=mlt_clim*rof_clim*snf_clim*rnf_clim
mltB=mlt_clim[ ~np.isnan(tmp) ]
rofB=rof_clim[ ~np.isnan(tmp) ]
snfB=snf_clim[ ~np.isnan(tmp) ]
rnfB=rnf_clim[ ~np.isnan(tmp) ]

#====================
# KDE fit :

xKDE=np.arange(0,130,0.1)
Xsigma=1.0

fit=np.zeros(np.shape(xKDE))

for kk in np.arange(0,np.size(xKDE)):
  kernel = np.exp(-((mltB+rnfB)/(snfB-mltB) - xKDE[kk]) ** 2 / (2 * Xsigma ** 2))
  kernel = kernel / sum(kernel)
  fit[kk] = sum( rofB * kernel )

#====================
# % of points experiencing runoff > 0.01 m/yr

thres=0.001

n1=0
n2=0
dx=1.0
xx=np.arange(0,130,dx)
Nrof=np.zeros((np.size(xx)))
Ndry=np.zeros((np.size(xx)))
for kk in np.arange(np.size(xx)):
  for pp in np.arange(np.size(rofB)):
     if (mltB[pp]+rnfB[pp])/(snfB[pp]-mltB[pp]) > xx[kk]-0.5*dx and (mltB[pp]+rnfB[pp])/(snfB[pp]-mltB[pp]) <= xx[kk]+0.5*dx:
       if rofB[pp] >= thres:
         Nrof[kk]=Nrof[kk]+1
       else:
         Ndry[kk]=Ndry[kk]+1
  if Ndry[kk]==0 and Nrof[kk]==0:
    Nrof[kk]=1
  if Nrof[kk]/(Nrof[kk]+Ndry[kk]) > 0.1 and n1==0:
    print xx[kk-1], Nrof[kk-1]/(Nrof[kk-1]+Ndry[kk-1])
    print xx[kk], Nrof[kk]/(Nrof[kk]+Ndry[kk])
    n1=1
  elif Nrof[kk]/(Nrof[kk]+Ndry[kk]) > 0.5 and n2==0:
    print xx[kk-1], Nrof[kk-1]/(Nrof[kk-1]+Ndry[kk-1])
    print xx[kk], Nrof[kk]/(Nrof[kk]+Ndry[kk])
    n2=1
print Nrof/(Nrof+Ndry)

#====================
fig, ax = plt.subplots()

print np.nanmin((mltB+rnfB)/(snfB-mltB))
print np.nanmax((mltB+rnfB)/(snfB-mltB))

#ax.scatter(mltA/snfA,rofA,s=5,c='cornflowerblue',label='Present')
ax.plot([0.0,130.0],[0.0,0.0],'-',color='k',linewidth=0.5)
ax.scatter((mltB+rnfB)/(snfB-mltB),rofB,s=5,c='red',alpha=0.2)
ax.plot([0.30, 0.30],[-1.0, 1.0],'--',linewidth=0.6,color='firebrick') # line showing the limit where more than 10% of the points are above the 1mm/yr threshold
ax.plot([0.85, 0.85],[-1.0, 1.0],'--',linewidth=0.6,color='firebrick') # line showing the limit where more than 50% of the points are above the 1mm/yr threshold
ax.plot(xKDE,fit,'-',color='firebrick',linewidth=0.8)
plt.xlim(0.0,130.0)
plt.ylim(-0.03,0.31)
#ax.legend()
ax.set_xlabel('(melt+rainfall)/(snowfall-melt)',size=12)
ax.set_ylabel('runoff (m.w.e/yr)',size=12)
ax.text(0.32,0.28,'10%',color='firebrick',fontsize=10)
ax.text(0.87,0.28,'50%',color='firebrick',fontsize=10)

#plt.show()
fig.savefig('runoff_vs_melt_on_snf_ratio_rain.pdf')
