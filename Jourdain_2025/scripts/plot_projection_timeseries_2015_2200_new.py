import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os

#NB: all files actually go to 2199, not 2200...

filez_grd='ensemble_proj_asmb_grd.npz'
filez_isf='ensemble_proj_asmb_isf.npz'

fig1, axs = plt.subplots(nrows=1,ncols=1,figsize=(18.0,11.0))
fig2, ays = plt.subplots(nrows=1,ncols=1,figsize=(18.0,11.0))
#axs = axs.ravel()

def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

grd=xr.open_dataset('RCM_ice_regrid_04000m.nc2',decode_cf=False)
msk_isf = ( grd.ICE_MAR - grd.GROUND ) * grd.af2
msk_grd = grd.GROUND * grd.af2

print('Total grounded ice sheet area = ',msk_grd.sum(dim=["x","y"]).values * 4. * 4. * 1.e-6,' million km2')
print('Total ice shelf area = ',msk_isf.sum(dim=["x","y"]).values * 4. * 4. * 1.e-6,' million km2')

# conversion to Gt/yr
fac=4.e3*4.e3*1.e-12*86400*365.25

model  = [ 'ACCESS-CM2', 'ACCESS-ESM1-5', 'CanESM5'   , 'CESM2-WACCM', 'GISS-E2-1-H', 'IPSL-CM6A-LR', 'MRI-ESM2-0', 'UKESM1-0-LL' ]
member = [ 'r1i1p1f1'  , 'r1i1p1f1'     , 'r1i1p1f1'  , 'r1i1p1f1'   , 'r1i1p1f2'   ,  'r1i1p1f1'   ,  'r1i1p1f1' , 'r4i1p1f2'    ]
color  = [ 'lightpink' , 'fuchsia'      , 'mediumblue', 'orangered'  , 'gold'       ,  'black'      , 'olivedrab' , 'darkgrey'    ] 

Nmod = np.size(model)

time_smoothed = moving_average(np.arange(1980,2200),21)
Ntime=2199-1980+1

############################################################

grd_ssp126 = np.zeros((Ntime,Nmod))
grd_ssp585 = np.zeros((Ntime,Nmod))
isf_ssp126 = np.zeros((Ntime,Nmod))
isf_ssp585 = np.zeros((Ntime,Nmod))
   
for kmod in np.arange(Nmod):
   
  print(model[kmod])
     
  file_his = 'MAR-'+model[kmod]+'_asmb_1980-2014_histo_regrid_04000m.nc'
  if not os.path.exists(file_his):
    file_his = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_1980-2014_histo_regrid_04000m_FROM_6_MODELS.nc'
  if ( model[kmod] == 'UKESM1-0-LL' ):
    file_his = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_1980-2014_histo_regrid_04000m_FROM_UKESM1-0-LL-r1i1p1f2-histo.nc'
   
  file_ssp126 = 'MAR-'+model[kmod]+'_asmb_2015-2200_ssp126_regrid_04000m.nc'
  if not os.path.exists(file_ssp126):
    file_ssp126 = 'MAR-'+model[kmod]+'_asmb_2015-2200_ssp126_regrid_04000m_FROM_ssp585.nc'
    if not os.path.exists(file_ssp126):
      file_ssp126 = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_2015-2200_ssp126_regrid_04000m_MERGED.nc'
     
  file_ssp585 = 'MAR-'+model[kmod]+'_asmb_2015-2200_ssp585_regrid_04000m.nc'
  if not os.path.exists(file_ssp585):
    file_ssp585 = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_2015-2200_ssp585_regrid_04000m_MERGED.nc'
     
  d1=xr.open_dataset(file_his,decode_cf=False)
  tmp = d1.asmb * msk_grd * fac ;
  txp = tmp.sum(dim=["x","y"]).values
  print(np.shape(txp),np.shape(grd_ssp126))
  grd_ssp126[0:2014-1980+1,kmod] = txp
  grd_ssp585[0:2014-1980+1,kmod] = txp
     
  d2=xr.open_dataset(file_ssp126,decode_cf=False)
  tmp = d2.asmb * msk_grd * fac ;
  txp = tmp.sum(dim=["x","y"]).values
  print(np.shape(txp),np.shape(grd_ssp126))
  grd_ssp126[2015-1980:2199-1980+1,kmod] = txp
     
  d3=xr.open_dataset(file_ssp585,decode_cf=False)
  tmp = d3.asmb * msk_grd * fac ;
  txp = tmp.sum(dim=["x","y"]).values
  grd_ssp585[2015-1980:2199-1980+1,kmod] = txp
  
  d4=xr.open_dataset(file_his,decode_cf=False)
  tmp = d4.asmb * msk_isf * fac ;
  txp = tmp.sum(dim=["x","y"]).values
  isf_ssp126[0:2014-1980+1,kmod] = txp
  isf_ssp585[0:2014-1980+1,kmod] = txp

  d5=xr.open_dataset(file_ssp126,decode_cf=False)
  tmp = d5.asmb * msk_isf * fac ;
  txp = tmp.sum(dim=["x","y"]).values
  isf_ssp126[2015-1980:2199-1980+1,kmod] = txp

  d6=xr.open_dataset(file_ssp585,decode_cf=False)
  tmp = d6.asmb * msk_isf * fac ;
  txp = tmp.sum(dim=["x","y"]).values
  isf_ssp585[2015-1980:2199-1980+1,kmod] = txp

  #==================================================================================

  axs.plot(time_smoothed,moving_average(grd_ssp126[:,kmod],21),color=color[kmod],linewidth=0.6,label=model[kmod]+' (SSP1-2.6)')
  axs.plot(time_smoothed,moving_average(grd_ssp585[:,kmod],21),color=color[kmod],linewidth=2.6,label=model[kmod]+' (SSP5-8.5)')

  ays.plot(time_smoothed,moving_average(isf_ssp126[:,kmod],21),color=color[kmod],linewidth=0.6,label=model[kmod]+' (SSP1-2.6)')
  ays.plot(time_smoothed,moving_average(isf_ssp585[:,kmod],21),color=color[kmod],linewidth=2.6,label=model[kmod]+' (SSP5-8.5)')

  # Sea level equivalent over 2000-2200 w.r.t. 1995-2014 :
  slr126 = ( np.sum(grd_ssp126[2000-1989:2200-1980,kmod])-np.sum(grd_ssp126[1995-1980:2015-1980,kmod]) ) / -361.8
  slr585 = ( np.sum(grd_ssp585[2000-1980:2200-1980,kmod])-np.sum(grd_ssp585[1995-1980:2015-1980,kmod]) ) / -361.8
  print('sea level rise (mm) ',model[kmod],slr126, slr585)

## Add weighted 16-model very likely range:
zg=np.load('ensemble_proj_asmb_grd.npz')
grd_pct05_ssp585 = zg['grd_pct05_ssp585'][1990-1860:2090-1860+1]
grd_pct95_ssp585 = zg['grd_pct95_ssp585'][1990-1860:2090-1860+1]
grd_pct05_ssp126 = zg['grd_pct05_ssp126'][1990-1860:2090-1860+1]
grd_pct95_ssp126 = zg['grd_pct95_ssp126'][1990-1860:2090-1860+1]
Nhalf = np.size(grd_pct05_ssp585)
axs.fill_between(time_smoothed[0:Nhalf],grd_pct05_ssp585,grd_pct95_ssp585,color='firebrick',alpha=0.2,label='weighted SSP5-8.5 range (5-95th pct)')
axs.fill_between(time_smoothed[0:Nhalf],grd_pct05_ssp126,grd_pct95_ssp126,color='cornflowerblue',alpha=0.3,label='weighted SSP1-2.6 range (5-95th pct)')

zi=np.load('ensemble_proj_asmb_isf.npz')
isf_pct05_ssp585 = zi['isf_pct05_ssp585'][1990-1860:2090-1860+1]
isf_pct95_ssp585 = zi['isf_pct95_ssp585'][1990-1860:2090-1860+1]
isf_pct05_ssp126 = zi['isf_pct05_ssp126'][1990-1860:2090-1860+1]
isf_pct95_ssp126 = zi['isf_pct95_ssp126'][1990-1860:2090-1860+1]
ays.fill_between(time_smoothed[0:Nhalf],isf_pct05_ssp585,isf_pct95_ssp585,color='firebrick',alpha=0.2,label='weighted SSP5-8.5 range (5-95th pct)')
ays.fill_between(time_smoothed[0:Nhalf],isf_pct05_ssp126,isf_pct95_ssp126,color='cornflowerblue',alpha=0.3,label='weighted SSP1-2.6 range (5-95th pct)')

## Add limit for negative SMB (MAR-RACMO-HIRHM range in Mottram et al., 2021):
axs.fill_between([1980,2200],[-2323,-2323],[-1929,-1929],hatch='\\/',facecolor='lightgrey')
ays.fill_between([1980,2200],[-471,-471],[-395,-395],hatch='\\/',facecolor='lightgrey',label='zero SMB limit')

#-
axs.set_title('SMB anomaly over the grounded ice sheet',fontsize=22,fontweight='bold')
axs.set_ylabel('Gt/yr',fontsize=20)
axs.set_xlim([2000,2200])
axs.tick_params(axis='x', labelsize=20)
axs.tick_params(axis='y', labelsize=20)
axs.grid(color = 'black', linestyle = 'dotted', linewidth=1)

ays.set_title('SMB anomaly over the ice shelves',fontsize=22,fontweight='bold')
ays.set_ylabel('Gt/yr',fontsize=20)
ays.set_xlim([2000,2200])
ays.tick_params(axis='x', labelsize=20)
ays.tick_params(axis='y', labelsize=20)
ays.grid(color = 'black', linestyle = 'dotted', linewidth=1)
ays.legend(loc='lower left',fontsize=16)

##########################################################################

fig1.savefig("projection_timeseries_8_models_2000_2200_grounded.pdf")
fig2.savefig("projection_timeseries_8_models_2000_2200_iceshelf.pdf")
