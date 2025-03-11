import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os

filez_grd='ensemble_proj_asmb_grd.npz'
filez_isf='ensemble_proj_asmb_isf.npz'

fig1, axs = plt.subplots(nrows=1,ncols=1,figsize=(18.0,11.0))
fig2, ays = plt.subplots(nrows=1,ncols=1,figsize=(18.0,11.0))
#axs = axs.ravel()

def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float,axis=0)
    ret[n:,:] = ret[n:,:] - ret[:-n,:]
    return ret[n - 1:,:] / n

def moving_average_1d(a, n=3):
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

model  = [ 'ACCESS-CM2', 'ACCESS-ESM1-5', 'CanESM5' , 'CESM2'    , 'CESM2-WACCM', 'CNRM-CM6-1', 'CNRM-ESM2-1', 'GFDL-CM4', 'GFDL-ESM4', 'GISS-E2-1-H', 'INM-CM5-0', 'IPSL-CM6A-LR', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'NorESM2-MM', 'UKESM1-0-LL' ]
member = [ 'r1i1p1f1'  , 'r1i1p1f1'     , 'r1i1p1f1', 'r11i1p1f1', 'r1i1p1f1'   , 'r1i1p1f2'  , 'r1i1p1f2'   , 'r1i1p1f1', 'r1i1p1f1' , 'r1i1p1f2'   , 'r1i1p1f1' , 'r1i1p1f1'    , 'r1i1p1f1'     , 'r1i1p1f1'  , 'r1i1p1f1'  , 'r1i1p1f2'    ]
weight = [      11     ,      24        ,       3   ,        6   ,      10      ,      10     ,      10      ,      24   ,      47    ,      41      ,      18    ,      12       ,      43        ,      39    ,       47     ,       5       ]
Nmod = np.size(model)
Nwgt = np.sum(weight)
Ntime = 2100-1850+1
print(Nwgt)

############################################################
if not os.path.exists(filez_grd):

   grd_ens_ssp126 = np.zeros((Ntime,Nwgt))
   grd_ens_ssp245 = np.zeros((Ntime,Nwgt))
   grd_ens_ssp585 = np.zeros((Ntime,Nwgt))
   
   Nwgtm1=0
   
   for kmod in np.arange(Nmod):
   
      print(model[kmod])
   
      file_1850 = 'MAR-'+model[kmod]+'_asmb_1850-1979_histo_regrid_04000m_EXTENDED.nc'
      if not os.path.exists(file_1850):
        file_1850 = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_1850-1979_histo_regrid_04000m_FROM_6_MODELS_EXTENDED.nc'
   
      file_his = 'MAR-'+model[kmod]+'_asmb_1980-2014_histo_regrid_04000m.nc'
      if not os.path.exists(file_his):
        file_his = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_1980-2014_histo_regrid_04000m_FROM_6_MODELS.nc'
   
      file_ssp126 = 'MAR-'+model[kmod]+'_asmb_2015-2100_ssp126_regrid_04000m.nc'
      if not os.path.exists(file_ssp126):
        file_ssp126 = 'MAR-'+model[kmod]+'_asmb_2015-2100_ssp126_regrid_04000m_FROM_ssp585.nc'
        if not os.path.exists(file_ssp126):
          file_ssp126 = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_2015-2100_ssp126_regrid_04000m_FROM_6_MODELS.nc'
   
      file_ssp245 = 'MAR-'+model[kmod]+'_asmb_2015-2100_ssp245_regrid_04000m.nc'
      if not os.path.exists(file_ssp245):
        file_ssp245 = 'MAR-'+model[kmod]+'_asmb_2015-2100_ssp245_regrid_04000m_FROM_ssp585.nc'
        if not os.path.exists(file_ssp245):
          file_ssp245 = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_2015-2100_ssp245_regrid_04000m_FROM_6_MODELS.nc'
   
      file_ssp585 = 'MAR-'+model[kmod]+'_asmb_2015-2100_ssp585_regrid_04000m.nc'
      if not os.path.exists(file_ssp585):
        file_ssp585 = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_2015-2100_ssp585_regrid_04000m_FROM_6_MODELS.nc'
   
      d1=xr.open_dataset(file_1850,decode_cf=False)
      tmp = d1.asmb * msk_grd * fac ;
      txp = tmp.sum(dim=["x","y"]).values
      for kwgt in np.arange(Nwgtm1,Nwgtm1+weight[kmod],1):
         grd_ens_ssp126[0:1979-1850+1,kwgt] = txp
         grd_ens_ssp245[0:1979-1850+1,kwgt] = txp
         grd_ens_ssp585[0:1979-1850+1,kwgt] = txp
   
      d2=xr.open_dataset(file_his,decode_cf=False)
      tmp = d2.asmb * msk_grd * fac ;
      txp = tmp.sum(dim=["x","y"]).values
      for kwgt in np.arange(Nwgtm1,Nwgtm1+weight[kmod],1):
         grd_ens_ssp126[1980-1850:2014-1850+1,kwgt] = txp
         grd_ens_ssp245[1980-1850:2014-1850+1,kwgt] = txp
         grd_ens_ssp585[1980-1850:2014-1850+1,kwgt] = txp
   
      if os.path.exists(file_ssp126):
        d3=xr.open_dataset(file_ssp126,decode_cf=False)
        tmp = d3.asmb * msk_grd * fac ;
        txp = tmp.sum(dim=["x","y"]).values
        for kwgt in np.arange(Nwgtm1,Nwgtm1+weight[kmod],1):
           grd_ens_ssp126[2015-1850:2100-1850+1,kwgt] = txp
      else: # GFDL-CM4 (no ssp126)
        for kwgt in np.arange(Nwgtm1,Nwgtm1+weight[kmod],1):
           grd_ens_ssp126[2015-1850:2100-1850+1,kwgt] = np.nan
   
      d4=xr.open_dataset(file_ssp245,decode_cf=False)
      tmp = d4.asmb * msk_grd * fac ;
      txp = tmp.sum(dim=["x","y"]).values
      for kwgt in np.arange(Nwgtm1,Nwgtm1+weight[kmod],1):
         grd_ens_ssp245[2015-1850:2100-1850+1,kwgt] = txp
   
      d5=xr.open_dataset(file_ssp585,decode_cf=False)
      tmp = d5.asmb * msk_grd * fac ;
      txp = tmp.sum(dim=["x","y"]).values
      for kwgt in np.arange(Nwgtm1,Nwgtm1+weight[kmod],1):
         grd_ens_ssp585[2015-1850:2100-1850+1,kwgt] = txp
   
      Nwgtm1=Nwgtm1+weight[kmod]
  
   txp = moving_average(grd_ens_ssp126,21)
   grd_pct05_ssp126 = np.nanpercentile(txp, 5.00,axis=1)
   grd_pct17_ssp126 = np.nanpercentile(txp,16.67,axis=1)
   grd_pct50_ssp126 = np.nanpercentile(txp,50.00,axis=1)
   grd_pct83_ssp126 = np.nanpercentile(txp,83.33,axis=1)
   grd_pct95_ssp126 = np.nanpercentile(txp,95.00,axis=1)
   grd_mean_ssp126  = np.nanmean(txp,axis=1)

   txp = moving_average(grd_ens_ssp245,21)
   grd_pct05_ssp245 = np.nanpercentile(txp, 5.00,axis=1)
   grd_pct17_ssp245 = np.nanpercentile(txp,16.67,axis=1)
   grd_pct50_ssp245 = np.nanpercentile(txp,50.00,axis=1)
   grd_pct83_ssp245 = np.nanpercentile(txp,83.33,axis=1)
   grd_pct95_ssp245 = np.nanpercentile(txp,95.00,axis=1)
   grd_mean_ssp245  = np.nanmean(txp,axis=1)

   txp = moving_average(grd_ens_ssp585,21)
   grd_pct05_ssp585 = np.nanpercentile(txp, 5.00,axis=1)
   grd_pct17_ssp585 = np.nanpercentile(txp,16.67,axis=1)
   grd_pct50_ssp585 = np.nanpercentile(txp,50.00,axis=1)
   grd_pct83_ssp585 = np.nanpercentile(txp,83.33,axis=1)
   grd_pct95_ssp585 = np.nanpercentile(txp,95.00,axis=1)
   grd_mean_ssp585  = np.nanmean(txp,axis=1)

   # 1900 to 2010 SLR (w.r.t. 1891-1910) 
   integral_1900_2010_mmSLR = - ( np.sum(grd_ens_ssp126[1900-1850:2010-1850+1,:],axis=0) - (2011-1900)*np.mean(grd_ens_ssp126[1891-1850:1910-1850+1,:],axis=0) ) / 361.8
   int_1900_2010_mmSLR_pct = np.nanpercentile(integral_1900_2010_mmSLR,[5.00,16.67,50.00,83.33,95.00])

   # 2000 to 2099 SLR (w.r.t. 1995-2014)
   integral_2000_2099_ssp126_mmSLR = - ( np.sum( grd_ens_ssp126[2000-1850:2099-1850+1,:],axis=0) - (2100-2000)*np.mean(grd_ens_ssp126[1995-1850:2014-1850+1,:],axis=0) ) / 361.8
   integral_2000_2099_ssp245_mmSLR = - ( np.sum( grd_ens_ssp245[2000-1850:2099-1850+1,:],axis=0) - (2100-2000)*np.mean(grd_ens_ssp245[1995-1850:2014-1850+1,:],axis=0) ) / 361.8
   integral_2000_2099_ssp585_mmSLR = - ( np.sum( grd_ens_ssp585[2000-1850:2099-1850+1,:],axis=0) - (2100-2000)*np.mean(grd_ens_ssp585[1995-1850:2014-1850+1,:],axis=0) ) / 361.8
   int_2000_2099_ssp126_mmSLR_pct = np.nanpercentile(integral_2000_2099_ssp126_mmSLR,[5.00,16.67,50.00,83.33,95.00])
   int_2000_2099_ssp245_mmSLR_pct = np.nanpercentile(integral_2000_2099_ssp245_mmSLR,[5.00,16.67,50.00,83.33,95.00])
   int_2000_2099_ssp585_mmSLR_pct = np.nanpercentile(integral_2000_2099_ssp585_mmSLR,[5.00,16.67,50.00,83.33,95.00])

   np.savez(filez_grd,\
   int_1900_2010_mmSLR_pct = int_1900_2010_mmSLR_pct, \
   int_2000_2099_ssp126_mmSLR_pct = int_2000_2099_ssp126_mmSLR_pct, \
   int_2000_2099_ssp245_mmSLR_pct = int_2000_2099_ssp245_mmSLR_pct, \
   int_2000_2099_ssp585_mmSLR_pct = int_2000_2099_ssp585_mmSLR_pct, \
   grd_pct05_ssp126 = grd_pct05_ssp126, \
   grd_pct17_ssp126 = grd_pct17_ssp126, \
   grd_pct50_ssp126 = grd_pct50_ssp126, \
   grd_pct83_ssp126 = grd_pct83_ssp126, \
   grd_pct95_ssp126 = grd_pct95_ssp126, \
   grd_mean_ssp126  = grd_mean_ssp126, \
   grd_pct05_ssp245 = grd_pct05_ssp245, \
   grd_pct17_ssp245 = grd_pct17_ssp245, \
   grd_pct50_ssp245 = grd_pct50_ssp245, \
   grd_pct83_ssp245 = grd_pct83_ssp245, \
   grd_pct95_ssp245 = grd_pct95_ssp245, \
   grd_mean_ssp245  = grd_mean_ssp245, \
   grd_pct05_ssp585 = grd_pct05_ssp585, \
   grd_pct17_ssp585 = grd_pct17_ssp585, \
   grd_pct50_ssp585 = grd_pct50_ssp585, \
   grd_pct83_ssp585 = grd_pct83_ssp585, \
   grd_pct95_ssp585 = grd_pct95_ssp585, \
   grd_mean_ssp585  = grd_mean_ssp585)

else:

   print('WARNING : starting from existing npz file : '+filez_grd)
   zz=np.load(filez_grd)
   int_1900_2010_mmSLR_pct = np.array(zz['int_1900_2010_mmSLR_pct'])
   int_2000_2099_ssp126_mmSLR_pct = np.array(zz['int_2000_2099_ssp126_mmSLR_pct'])
   int_2000_2099_ssp245_mmSLR_pct = np.array(zz['int_2000_2099_ssp245_mmSLR_pct'])
   int_2000_2099_ssp585_mmSLR_pct = np.array(zz['int_2000_2099_ssp585_mmSLR_pct'])
   grd_pct05_ssp126 = zz['grd_pct05_ssp126']
   grd_pct17_ssp126 = zz['grd_pct17_ssp126']
   grd_pct50_ssp126 = zz['grd_pct50_ssp126']
   grd_pct83_ssp126 = zz['grd_pct83_ssp126']
   grd_pct95_ssp126 = zz['grd_pct95_ssp126']
   grd_mean_ssp126  = zz['grd_mean_ssp126']
   grd_pct05_ssp245 = zz['grd_pct05_ssp245']
   grd_pct17_ssp245 = zz['grd_pct17_ssp245']
   grd_pct50_ssp245 = zz['grd_pct50_ssp245']
   grd_pct83_ssp245 = zz['grd_pct83_ssp245']
   grd_pct95_ssp245 = zz['grd_pct95_ssp245']
   grd_mean_ssp245  = zz['grd_mean_ssp245']
   grd_pct05_ssp585 = zz['grd_pct05_ssp585']
   grd_pct17_ssp585 = zz['grd_pct17_ssp585']
   grd_pct50_ssp585 = zz['grd_pct50_ssp585']
   grd_pct83_ssp585 = zz['grd_pct83_ssp585']
   grd_pct95_ssp585 = zz['grd_pct95_ssp585']
   grd_mean_ssp585  = zz['grd_mean_ssp585']

#==================================================================================

print('SLR 1900-2010 = ',round(int_1900_2010_mmSLR_pct[2],2),'   (',round(int_1900_2010_mmSLR_pct[1],2),' -- ',round(int_1900_2010_mmSLR_pct[3],2),')   [',round(int_1900_2010_mmSLR_pct[0],2),' -- ',round(int_1900_2010_mmSLR_pct[4],2),']   mm SLR')
print('SLR 2000-2099 (SSP126) = ',round(int_2000_2099_ssp126_mmSLR_pct[2],2),'   (',round(int_2000_2099_ssp126_mmSLR_pct[1],2),' -- ',round(int_2000_2099_ssp126_mmSLR_pct[3],2),')   [',round(int_2000_2099_ssp126_mmSLR_pct[0],2),' -- ',round(int_2000_2099_ssp126_mmSLR_pct[4],2),']   mm SLR')
print('SLR 2000-2099 (SSP245) = ',round(int_2000_2099_ssp245_mmSLR_pct[2],2),'   (',round(int_2000_2099_ssp245_mmSLR_pct[1],2),' -- ',round(int_2000_2099_ssp245_mmSLR_pct[3],2),')   [',round(int_2000_2099_ssp245_mmSLR_pct[0],2),' -- ',round(int_2000_2099_ssp245_mmSLR_pct[4],2),']   mm SLR')
print('SLR 2000-2099 (SSP585) = ',round(int_2000_2099_ssp585_mmSLR_pct[2],2),'   (',round(int_2000_2099_ssp585_mmSLR_pct[1],2),' -- ',round(int_2000_2099_ssp585_mmSLR_pct[3],2),')   [',round(int_2000_2099_ssp585_mmSLR_pct[0],2),' -- ',round(int_2000_2099_ssp585_mmSLR_pct[4],2),']   mm SLR')

time_smoothed = moving_average_1d(np.arange(1850,2101),21)

axs.fill_between(time_smoothed,grd_pct05_ssp585,grd_pct95_ssp585,color='firebrick',alpha=0.1,label='SSP5-8.5 very likely range (5-95th pct)')
axs.plot(time_smoothed,grd_pct50_ssp585,color='firebrick',linewidth=3,label='SSP5-8.5 median')

axs.fill_between(time_smoothed,grd_pct05_ssp245,grd_pct95_ssp245,color='orange',alpha=0.1,label='SSP2-4.5 very likely range (5-95th pct)')
axs.plot(time_smoothed,grd_pct50_ssp245,color='orange',linewidth=3,label='SSP2-4.5 median')

axs.fill_between(time_smoothed,grd_pct05_ssp126,grd_pct95_ssp126,color='cornflowerblue',alpha=0.1,label='SSP1-2.6 very likely range (5-95th pct)')
axs.plot(time_smoothed,grd_pct50_ssp126,color='cornflowerblue',linewidth=3,label='SSP1-2.6 median')

axs.plot(time_smoothed[0:2014-1850-20],grd_pct50_ssp126[0:2014-1850-20],color='gray',linewidth=4,label='historical median')
# 66% confidence interval as error bars:
st=[0, 40, 90, 140]
xt=time_smoothed[st]; print(xt)
ym=grd_pct50_ssp126[st] ; print(ym)
yl=ym-grd_pct17_ssp126[st]; print(yl)
yu=grd_pct83_ssp126[st]-ym; print(yu)
er=[yl,yu]
axs.errorbar(xt,ym,yerr=er,fmt='o',color='gray',elinewidth=2,markersize=12, capsize=8,label='likely range (17-83th pct)')
#-
xtb=time_smoothed[-1]; print(xtb)
ym126=grd_pct50_ssp126[-1]; print('ym126 = ',ym126)
yl126=ym126-grd_pct17_ssp126[-1]
yu126=grd_pct83_ssp126[-1]-ym126
ym245=grd_pct50_ssp245[-1]; print('ym245= ',ym245)
yl245=ym245-grd_pct17_ssp245[-1]
yu245=grd_pct83_ssp245[-1]-ym245
ym585=grd_pct50_ssp585[-1]; print('ym585= ',ym585)
yl585=ym585-grd_pct17_ssp585[-1]
yu585=grd_pct83_ssp585[-1]-ym585
axs.errorbar(xtb,ym126,yerr=[[yl126],[yu126]],fmt='o',color='cornflowerblue',elinewidth=2,markersize=12, capsize=8)
axs.errorbar(xtb-0.75,ym245,yerr=[[yl245],[yu245]],fmt='o',color='orange',elinewidth=2,markersize=12, capsize=8)
axs.errorbar(xtb,ym585,yerr=[[yl585],[yu585]],fmt='o',color='firebrick',elinewidth=2,markersize=12, capsize=8)

axs.set_title('SMB anomaly over the grounded ice sheet',fontsize=22,fontweight='bold')
axs.set_ylabel('Gt/yr',fontsize=20)
axs.legend(loc='upper left',fontsize=20)
axs.set_xlim([1850,2100])
axs.set_ylim([-220,750])
axs.tick_params(axis='x', labelsize=20)
axs.tick_params(axis='y', labelsize=20)
axs.grid(color = 'black', linestyle = 'dotted', linewidth=1)

############################################################
if not os.path.exists(filez_isf):

   isf_ens_ssp126 = np.zeros((Ntime,Nwgt))
   isf_ens_ssp245 = np.zeros((Ntime,Nwgt))
   isf_ens_ssp585 = np.zeros((Ntime,Nwgt))
   
   Nwgtm1=0
   
   for kmod in np.arange(Nmod):
   
      print(model[kmod])
   
      file_1850 = 'MAR-'+model[kmod]+'_asmb_1850-1979_histo_regrid_04000m_EXTENDED.nc'
      if not os.path.exists(file_1850):
        file_1850 = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_1850-1979_histo_regrid_04000m_FROM_6_MODELS_EXTENDED.nc'
   
      file_his = 'MAR-'+model[kmod]+'_asmb_1980-2014_histo_regrid_04000m.nc'
      if not os.path.exists(file_his):
        file_his = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_1980-2014_histo_regrid_04000m_FROM_6_MODELS.nc'
   
      file_ssp126 = 'MAR-'+model[kmod]+'_asmb_2015-2100_ssp126_regrid_04000m.nc'
      if not os.path.exists(file_ssp126):
        file_ssp126 = 'MAR-'+model[kmod]+'_asmb_2015-2100_ssp126_regrid_04000m_FROM_ssp585.nc'
        if not os.path.exists(file_ssp126):
          file_ssp126 = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_2015-2100_ssp126_regrid_04000m_FROM_6_MODELS.nc'
   
      file_ssp245 = 'MAR-'+model[kmod]+'_asmb_2015-2100_ssp245_regrid_04000m.nc'
      if not os.path.exists(file_ssp245):
        file_ssp245 = 'MAR-'+model[kmod]+'_asmb_2015-2100_ssp245_regrid_04000m_FROM_ssp585.nc'
        if not os.path.exists(file_ssp245):
          file_ssp245 = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_2015-2100_ssp245_regrid_04000m_FROM_6_MODELS.nc'
   
      file_ssp585 = 'MAR-'+model[kmod]+'_asmb_2015-2100_ssp585_regrid_04000m.nc'
      if not os.path.exists(file_ssp585):
        file_ssp585 = 'MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_2015-2100_ssp585_regrid_04000m_FROM_6_MODELS.nc'
   
      d1=xr.open_dataset(file_1850,decode_cf=False)
      tmp = d1.asmb * msk_isf * fac ;
      txp = tmp.sum(dim=["x","y"]).values
      for kwgt in np.arange(Nwgtm1,Nwgtm1+weight[kmod],1):
         isf_ens_ssp126[0:1979-1850+1,kwgt] = txp
         isf_ens_ssp245[0:1979-1850+1,kwgt] = txp
         isf_ens_ssp585[0:1979-1850+1,kwgt] = txp
   
      d2=xr.open_dataset(file_his,decode_cf=False)
      tmp = d2.asmb * msk_isf * fac ;
      txp = tmp.sum(dim=["x","y"]).values
      for kwgt in np.arange(Nwgtm1,Nwgtm1+weight[kmod],1):
         isf_ens_ssp126[1980-1850:2014-1850+1,kwgt] = txp
         isf_ens_ssp245[1980-1850:2014-1850+1,kwgt] = txp
         isf_ens_ssp585[1980-1850:2014-1850+1,kwgt] = txp
   
      if os.path.exists(file_ssp126):
        d3=xr.open_dataset(file_ssp126,decode_cf=False)
        tmp = d3.asmb * msk_isf * fac ;
        txp = tmp.sum(dim=["x","y"]).values
        for kwgt in np.arange(Nwgtm1,Nwgtm1+weight[kmod],1):
           isf_ens_ssp126[2015-1850:2100-1850+1,kwgt] = txp
      else: # GFDL-CM4 (no ssp126)
        for kwgt in np.arange(Nwgtm1,Nwgtm1+weight[kmod],1):
           isf_ens_ssp126[2015-1850:2100-1850+1,kwgt] = np.nan
   
      d4=xr.open_dataset(file_ssp245,decode_cf=False)
      tmp = d4.asmb * msk_isf * fac ;
      txp = tmp.sum(dim=["x","y"]).values
      for kwgt in np.arange(Nwgtm1,Nwgtm1+weight[kmod],1):
         isf_ens_ssp245[2015-1850:2100-1850+1,kwgt] = txp
   
      d5=xr.open_dataset(file_ssp585,decode_cf=False)
      tmp = d5.asmb * msk_isf * fac ;
      txp = tmp.sum(dim=["x","y"]).values
      for kwgt in np.arange(Nwgtm1,Nwgtm1+weight[kmod],1):
         isf_ens_ssp585[2015-1850:2100-1850+1,kwgt] = txp
   
      Nwgtm1=Nwgtm1+weight[kmod]
   
   txp = moving_average(isf_ens_ssp126,21)
   isf_pct05_ssp126 = np.nanpercentile(txp, 5.00,axis=1)
   isf_pct17_ssp126 = np.nanpercentile(txp,16.67,axis=1)
   isf_pct50_ssp126 = np.nanpercentile(txp,50.00,axis=1)
   isf_pct83_ssp126 = np.nanpercentile(txp,83.33,axis=1)
   isf_pct95_ssp126 = np.nanpercentile(txp,95.00,axis=1)
   isf_mean_ssp126  = np.nanmean(txp,axis=1)
   
   txp = moving_average(isf_ens_ssp245,21)
   isf_pct05_ssp245 = np.nanpercentile(txp, 5.00,axis=1)
   isf_pct17_ssp245 = np.nanpercentile(txp,16.67,axis=1)
   isf_pct50_ssp245 = np.nanpercentile(txp,50.00,axis=1)
   isf_pct83_ssp245 = np.nanpercentile(txp,83.33,axis=1)
   isf_pct95_ssp245 = np.nanpercentile(txp,95.00,axis=1)
   isf_mean_ssp245  = np.nanmean(txp,axis=1)
  
   txp = moving_average(isf_ens_ssp585,21)
   isf_pct05_ssp585 = np.nanpercentile(txp, 5.00,axis=1)
   isf_pct17_ssp585 = np.nanpercentile(txp,16.67,axis=1)
   isf_pct50_ssp585 = np.nanpercentile(txp,50.00,axis=1)
   isf_pct83_ssp585 = np.nanpercentile(txp,83.33,axis=1)
   isf_pct95_ssp585 = np.nanpercentile(txp,95.00,axis=1)
   isf_mean_ssp585  = np.nanmean(txp,axis=1)
 
   np.savez(filez_isf,\
   isf_pct05_ssp126 = isf_pct05_ssp126, \
   isf_pct17_ssp126 = isf_pct17_ssp126, \
   isf_pct50_ssp126 = isf_pct50_ssp126, \
   isf_pct83_ssp126 = isf_pct83_ssp126, \
   isf_pct95_ssp126 = isf_pct95_ssp126, \
   isf_mean_ssp126  = isf_mean_ssp126, \
   isf_pct05_ssp245 = isf_pct05_ssp245, \
   isf_pct17_ssp245 = isf_pct17_ssp245, \
   isf_pct50_ssp245 = isf_pct50_ssp245, \
   isf_pct83_ssp245 = isf_pct83_ssp245, \
   isf_pct95_ssp245 = isf_pct95_ssp245, \
   isf_mean_ssp245  = isf_mean_ssp245, \
   isf_pct05_ssp585 = isf_pct05_ssp585, \
   isf_pct17_ssp585 = isf_pct17_ssp585, \
   isf_pct50_ssp585 = isf_pct50_ssp585, \
   isf_pct83_ssp585 = isf_pct83_ssp585, \
   isf_pct95_ssp585 = isf_pct95_ssp585, \
   isf_mean_ssp585  = isf_mean_ssp585)

else:

   print('WARNING : starting from existing npz file : '+filez_isf)
   zz=np.load(filez_isf)
   isf_pct05_ssp126 = zz['isf_pct05_ssp126']
   isf_pct17_ssp126 = zz['isf_pct17_ssp126']
   isf_pct50_ssp126 = zz['isf_pct50_ssp126']
   isf_pct83_ssp126 = zz['isf_pct83_ssp126']
   isf_pct95_ssp126 = zz['isf_pct95_ssp126']
   isf_mean_ssp126  = zz['isf_mean_ssp126']
   isf_pct05_ssp245 = zz['isf_pct05_ssp245']
   isf_pct17_ssp245 = zz['isf_pct17_ssp245']
   isf_pct50_ssp245 = zz['isf_pct50_ssp245']
   isf_pct83_ssp245 = zz['isf_pct83_ssp245']
   isf_pct95_ssp245 = zz['isf_pct95_ssp245']
   isf_mean_ssp245  = zz['isf_mean_ssp245']
   isf_pct05_ssp585 = zz['isf_pct05_ssp585']
   isf_pct17_ssp585 = zz['isf_pct17_ssp585']
   isf_pct50_ssp585 = zz['isf_pct50_ssp585']
   isf_pct83_ssp585 = zz['isf_pct83_ssp585']
   isf_pct95_ssp585 = zz['isf_pct95_ssp585']
   isf_mean_ssp585  = zz['isf_mean_ssp585']

#==================================================================================

time_smoothed = moving_average_1d(np.arange(1850,2101),21)

ays.fill_between(time_smoothed,isf_pct05_ssp585,isf_pct95_ssp585,color='firebrick',alpha=0.1,label='SSP5-8.5 very likely range (5-95th pct)')
ays.plot(time_smoothed,isf_pct50_ssp585,color='firebrick',linewidth=3,label='SSP5-8.5 median')

ays.fill_between(time_smoothed,isf_pct05_ssp245,isf_pct95_ssp245,color='orange',alpha=0.1,label='SSP2-4.5 very likely range (5-95th pct)')
ays.plot(time_smoothed,isf_pct50_ssp245,color='orange',linewidth=3,label='SSP2-4.5 median')

ays.fill_between(time_smoothed,isf_pct05_ssp126,isf_pct95_ssp126,color='cornflowerblue',alpha=0.1,label='SSP1-2.6 very likely range (5-95th pct)')
ays.plot(time_smoothed,isf_pct50_ssp126,color='cornflowerblue',linewidth=3,label='SSP1-2.6 median')

ays.plot(time_smoothed[0:2014-1850-20],isf_pct50_ssp126[0:2014-1850-20],color='gray',linewidth=4,label='historical median')
# 66% confidence interval as error bars:
st=[0, 40, 90, 140]
xt=time_smoothed[st]; print(xt)
ym=isf_pct50_ssp126[st]; print(ym)
yl=ym-isf_pct17_ssp126[st]
yu=isf_pct83_ssp126[st]-ym
er=[yl,yu]
ays.errorbar(xt,ym,yerr=er,fmt='o',color='gray',elinewidth=2,markersize=12, capsize=8,label='likely range (66.7%)')
#-
xtb=time_smoothed[-1]
ym126=isf_pct50_ssp126[-1]
yl126=ym126-isf_pct17_ssp126[-1]
yu126=isf_pct83_ssp126[-1]-ym126
ym245=isf_pct50_ssp245[-1]
yl245=ym245-isf_pct17_ssp245[-1]
yu245=isf_pct83_ssp245[-1]-ym245
ym585=isf_pct50_ssp585[-1]
yl585=ym585-isf_pct17_ssp585[-1]
yu585=isf_pct83_ssp585[-1]-ym585
ays.errorbar(xtb,ym126,yerr=[[yl126],[yu126]],fmt='o',color='cornflowerblue',elinewidth=2,markersize=12, capsize=8)
ays.errorbar(xtb-0.75,ym245,yerr=[[yl245],[yu245]],fmt='o',color='orange',elinewidth=2,markersize=12, capsize=8)
ays.errorbar(xtb,ym585,yerr=[[yl585],[yu585]],fmt='o',color='firebrick',elinewidth=2,markersize=12, capsize=8)

ays.set_title('SMB anomaly over the ice shelves',fontsize=22,fontweight='bold')
ays.set_ylabel('Gt/yr',fontsize=20)
ays.set_xlim([1850,2100])
ays.set_ylim([-380,35])
ays.tick_params(axis='x', labelsize=20)
ays.tick_params(axis='y', labelsize=20)
ays.grid(color = 'black', linestyle = 'dotted', linewidth=1)

##########################################################################

fig1.savefig("projection_timeseries_16_models_1850_2100_grounded.pdf")
fig2.savefig("projection_timeseries_16_models_1850_2100_iceshelf.pdf")
