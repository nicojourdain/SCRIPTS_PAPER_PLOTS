import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os
from scipy.stats import norm

timewin = 11 # nb of years over which runoff is averaged to be compared to threshold
#threshold = 100.0 # threshold over which hydrofracturing is possible (kg m-2 yr-1)

# threshold over which hydrofracturing is possible (kg m-2 yr-1)
# thresholds taken at 5%, 10%, 15%, 20%, 25%,... of a normal distribution centered on 150.0 with sigma=61.0
x=np.linspace(0,500,5001)
rv=norm(loc=150,scale=61) # distrib that gives 90% of the values between 50 and 250 kg m-2 yr-1 (not <50 to avoid too many favorable conditions in preind and 250 between 200 and 300).
cdf=rv.cdf(x)
per = np.arange(5,100,5)
Nwt = np.size(per)
threshold = np.zeros((Nwt))
for pct in np.arange(Nwt):
   for kk in np.arange(np.size(x)):
      if ( cdf[kk] >= per[pct]*0.01 ):
         threshold[pct]=x[kk]
         break
print('thresholds : ',threshold)

filez_runoff_dates='runoff_dates_2200_weight.npz'

fig, axs = plt.subplots(nrows=2,ncols=1,figsize=(18.0,22.0))
axs = axs.ravel()

# moving average along axis=0 for 3d array
def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float, axis=0)
    ret[n:,:,:] = ret[n:,:,:] - ret[:-n,:,:]
    return ret[n - 1:,:,:] / n

# moving average for 1d array
def moving_average_1d(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

# MAR masks:
grd=xr.open_dataset('RCM_ice_regrid_04000m.nc2',decode_cf=False)
msk_isf = ( grd.ICE_MAR - grd.GROUND ) * grd.af2
msk_grd = grd.GROUND * grd.af2

# Basin masks:
basin=xr.open_dataset('/data/njourdain/DATA_ISMIP6/Mask_Iceshelf_IMBIE_4km.nc')
Nbasin = basin.ID.size

# Fraction of the coastline occupied by ice shelves in each basin:
zratio=np.load('fraction_coast_isf.npz')
ratio=np.array(zratio['ratio'])

print('Total grounded ice sheet area = ',msk_grd.sum(dim=["x","y"]).values * 4. * 4. * 1.e-6,' million km2')
print('Total ice shelf area = ',msk_isf.sum(dim=["x","y"]).values * 4. * 4. * 1.e-6,' million km2')

model  = [ 'ACCESS-CM2', 'ACCESS-ESM1-5', 'CanESM5' , 'CESM2-WACCM', 'GISS-E2-1-H', 'IPSL-CM6A-LR', 'MRI-ESM2-0', 'UKESM1-0-LL' ]
member = [ 'r1i1p1f1'  , 'r1i1p1f1'     , 'r1i1p1f1', 'r1i1p1f1'   , 'r1i1p1f2'   , 'r1i1p1f1'    , 'r1i1p1f1'  , 'r4i1p1f2'    ]
weight = [      11     ,      24        ,       3   ,      10      ,      41      ,      12       ,       39    ,       5       ]
Nmod = np.size(model)
Nwgt = np.sum(weight)
Ntime = 2199-2015+1
year = moving_average_1d(np.arange(2015,2200,1).astype('int'),timewin)
print(Nwgt)

############################################################
if not os.path.exists(filez_runoff_dates):

   avgru_isf_ens_ssp126 = np.zeros((Ntime,Nwgt*Nwt,Nbasin)) * np.nan
   avgru_isf_ens_ssp245 = np.zeros((Ntime,Nwgt*Nwt,Nbasin)) * np.nan
   avgru_isf_ens_ssp585 = np.zeros((Ntime,Nwgt*Nwt,Nbasin)) * np.nan

   avgru_bas_ens_ssp126 = np.zeros((Ntime,Nwgt*Nwt,Nbasin)) * np.nan
   avgru_bas_ens_ssp245 = np.zeros((Ntime,Nwgt*Nwt,Nbasin)) * np.nan
   avgru_bas_ens_ssp585 = np.zeros((Ntime,Nwgt*Nwt,Nbasin)) * np.nan

   Nwgtm1=0
   
   for kmod in np.arange(Nmod):
   
      print(model[kmod])

      file_clim = 'MAR-'+model[kmod]+'_ru_1995-2014_clim_regrid_04000m.nc'
      if not os.path.exists(file_clim):
        file_clim = 'MAR-'+model[kmod]+'-'+member[kmod]+'_ru_1995-2014_clim_regrid_04000m_FROM_6_MODELS.nc' 
      print(file_clim)

      file_ssp126 = 'MAR-'+model[kmod]+'_aru_2015-2200_ssp126_regrid_04000m.nc'
      if not os.path.exists(file_ssp126):
        file_ssp126 = 'MAR-'+model[kmod]+'_aru_2015-2200_ssp126_regrid_04000m_FROM_ssp585.nc'
        if not os.path.exists(file_ssp126):
          file_ssp126 = 'MAR-'+model[kmod]+'-'+member[kmod]+'_aru_2015-2200_ssp126_regrid_04000m_MERGED.nc'
      print(file_ssp126)
        
      file_ssp585 = 'MAR-'+model[kmod]+'_aru_2015-2200_ssp585_regrid_04000m.nc'
      if not os.path.exists(file_ssp585):
        file_ssp585 = 'MAR-'+model[kmod]+'-'+member[kmod]+'_aru_2015-2200_ssp585_regrid_04000m_MERGED.nc'
      print(file_ssp585)

      d0=xr.open_dataset(file_clim,decode_cf=False)
      d3=xr.open_dataset(file_ssp126,decode_cf=False)
      d5=xr.open_dataset(file_ssp585,decode_cf=False)

      for kbasin in np.arange(Nbasin):
         mskisf = ( grd.ICE_MAR - grd.GROUND ) * ( basin.Iceshelf_extrap.where( (basin.Iceshelf_extrap == basin.ID[kbasin]) ) * 0 + 1) * grd.af2
         mskgrd = grd.GROUND * ( basin.Iceshelf_extrap.where( (basin.Iceshelf_extrap == basin.ID.values[kbasin]) ) * 0 + 1) * grd.af2
         mskbas = grd.AIS * ( basin.Iceshelf_extrap.where( (basin.Iceshelf_extrap == basin.ID[kbasin]) ) * 0 + 1) * grd.af2
         area_isf = mskisf.sum(skipna=True).values # in nb of grid points
         area_grd = mskgrd.sum(skipna=True).values
         if (   ( ( ( area_isf*4.0**2 > 1500 ) & ( area_grd*4.0**2 > 15000 ) ) \
              | ( basin.NAME.values[kbasin] == 'LarsenA' )   \
              | ( basin.NAME.values[kbasin] == 'LarsenB' )   \
              | ( basin.NAME.values[kbasin] == 'Wordie'  )   \
              | ( basin.NAME.values[kbasin] == 'Venable' )   \
              | ( basin.NAME.values[kbasin] == 'Cosgrove') ) \
              & ( not basin.NAME.values[kbasin] == '' ) ):
           # SSP126:
           tmp_isf = ((d0.ru+d3.aru) * mskisf).sum(dim=["x","y"]) * 365.25 * 86400 / area_isf # only runoff produced on the ice shelf [kg m-2 yr-1]
           tmp_bas = ((d0.ru+d3.aru) * ( ratio[kbasin]*mskbas + (1-ratio[kbasin])*mskisf ) ).sum(dim=["x","y"]) * 365.25 * 86400 / area_isf
           for kwgt in np.arange(Nwgtm1,Nwgtm1+weight[kmod]*Nwt,1):
              kkw = np.mod((kwgt-Nwgtm1),Nwt)
              avgru_isf_ens_ssp126[:,kwgt,kbasin] = tmp_isf - threshold[kkw]
              avgru_bas_ens_ssp126[:,kwgt,kbasin] = tmp_bas - threshold[kkw]
           # SSP585:
           tmp_isf = ((d0.ru+d5.aru) * mskisf).sum(dim=["x","y"]) * 365.25 * 86400 / area_isf # only runoff produced on the ice shelf [kg m-2 yr-1]
           tmp_bas = ((d0.ru+d5.aru) * ( ratio[kbasin]*mskbas + (1-ratio[kbasin])*mskisf ) ).sum(dim=["x","y"]) * 365.25 * 86400 / area_isf
           for kwgt in np.arange(Nwgtm1,Nwgtm1+weight[kmod]*Nwt,1):
              kkw = np.mod((kwgt-Nwgtm1),Nwt)
              avgru_isf_ens_ssp585[:,kwgt,kbasin] = tmp_isf - threshold[kkw]
              avgru_bas_ens_ssp585[:,kwgt,kbasin] = tmp_bas - threshold[kkw]
         else:
           for kwgt in np.arange(Nwgtm1,Nwgtm1+weight[kmod]*Nwt,1):
              kkw = np.mod((kwgt-Nwgtm1),Nwt)
              avgru_isf_ens_ssp126[:,kwgt,kbasin] = np.nan
              avgru_bas_ens_ssp126[:,kwgt,kbasin] = np.nan
              avgru_isf_ens_ssp585[:,kwgt,kbasin] = np.nan
              avgru_bas_ens_ssp585[:,kwgt,kbasin] = np.nan

      Nwgtm1=Nwgtm1+weight[kmod]*Nwt
  
   #-- Smoothing and weighted percentiles across models: 
   avgru_isf_ssp126_pct05 = np.nanpercentile(moving_average(avgru_isf_ens_ssp126, timewin), 5.00,axis=1)
   avgru_isf_ssp126_pct17 = np.nanpercentile(moving_average(avgru_isf_ens_ssp126, timewin),16.67,axis=1)
   avgru_isf_ssp126_pct83 = np.nanpercentile(moving_average(avgru_isf_ens_ssp126, timewin),83.33,axis=1)
   avgru_isf_ssp126_pct95 = np.nanpercentile(moving_average(avgru_isf_ens_ssp126, timewin),95.00,axis=1)

   avgru_isf_ssp585_pct05 = np.nanpercentile(moving_average(avgru_isf_ens_ssp585, timewin), 5.00,axis=1)
   avgru_isf_ssp585_pct17 = np.nanpercentile(moving_average(avgru_isf_ens_ssp585, timewin),16.67,axis=1)
   avgru_isf_ssp585_pct83 = np.nanpercentile(moving_average(avgru_isf_ens_ssp585, timewin),83.33,axis=1)
   avgru_isf_ssp585_pct95 = np.nanpercentile(moving_average(avgru_isf_ens_ssp585, timewin),95.00,axis=1)

   avgru_bas_ssp126_pct05 = np.nanpercentile(moving_average(avgru_bas_ens_ssp126, timewin), 5.00,axis=1)
   avgru_bas_ssp126_pct17 = np.nanpercentile(moving_average(avgru_bas_ens_ssp126, timewin),16.67,axis=1)
   avgru_bas_ssp126_pct83 = np.nanpercentile(moving_average(avgru_bas_ens_ssp126, timewin),83.33,axis=1)
   avgru_bas_ssp126_pct95 = np.nanpercentile(moving_average(avgru_bas_ens_ssp126, timewin),95.00,axis=1)

   avgru_bas_ssp585_pct05 = np.nanpercentile(moving_average(avgru_bas_ens_ssp585, timewin), 5.00,axis=1)
   avgru_bas_ssp585_pct17 = np.nanpercentile(moving_average(avgru_bas_ens_ssp585, timewin),16.67,axis=1)
   avgru_bas_ssp585_pct83 = np.nanpercentile(moving_average(avgru_bas_ens_ssp585, timewin),83.33,axis=1)
   avgru_bas_ssp585_pct95 = np.nanpercentile(moving_average(avgru_bas_ens_ssp585, timewin),95.00,axis=1)

   #-- Years from which the runoff threshold is crossed
   year_isf_ssp126_pct05 = np.zeros((Nbasin)).astype('int'); year_isf_ssp585_pct05 = np.zeros((Nbasin)).astype('int')
   year_isf_ssp126_pct17 = np.zeros((Nbasin)).astype('int'); year_isf_ssp585_pct17 = np.zeros((Nbasin)).astype('int')
   year_isf_ssp126_pct83 = np.zeros((Nbasin)).astype('int'); year_isf_ssp585_pct83 = np.zeros((Nbasin)).astype('int')
   year_isf_ssp126_pct95 = np.zeros((Nbasin)).astype('int'); year_isf_ssp585_pct95 = np.zeros((Nbasin)).astype('int')

   year_bas_ssp126_pct05 = np.zeros((Nbasin)).astype('int'); year_bas_ssp585_pct05 = np.zeros((Nbasin)).astype('int')
   year_bas_ssp126_pct17 = np.zeros((Nbasin)).astype('int'); year_bas_ssp585_pct17 = np.zeros((Nbasin)).astype('int')
   year_bas_ssp126_pct83 = np.zeros((Nbasin)).astype('int'); year_bas_ssp585_pct83 = np.zeros((Nbasin)).astype('int')
   year_bas_ssp126_pct95 = np.zeros((Nbasin)).astype('int'); year_bas_ssp585_pct95 = np.zeros((Nbasin)).astype('int')

   for kbasin in np.arange(Nbasin):

      for typ in ['isf', 'bas']:
         for scenar in ['ssp126', 'ssp585']:
            for pct in ['05', '17', '83', '95']:
               if ( eval("np.max(avgru_"+typ+"_"+scenar+"_pct"+pct+"[:,kbasin])") < 0. ):
                 # always under threshold -> set date to 2500
                 exec("year_"+typ+"_"+scenar+"_pct"+pct+"[kbasin] = 2500")
               elif (   ( eval("np.min(avgru_"+typ+"_"+scenar+"_pct"+pct+"[:,kbasin])") >= 0. ) \
                      | ( eval("avgru_"+typ+"_"+scenar+"_pct"+pct+"[0,kbasin]") >= 0. ) ):
                 # always above threshold or first point already above -> set date to 1500
                 exec("year_"+typ+"_"+scenar+"_pct"+pct+"[kbasin] = 1500")
               elif ( eval("np.sum(~np.isnan(avgru_"+typ+"_"+scenar+"_pct"+pct+"[:,kbasin]))") == 0 ):
                 # non calculated basin -> set to zero
                 exec("year_"+typ+"_"+scenar+"_pct"+pct+"[kbasin] = 0")
               else:
                 # date = end of time window:
                 dfp = eval("np.argwhere(avgru_"+typ+"_"+scenar+"_pct"+pct+"[:,kbasin] >= 0.)")
                 if ( np.size(dfp) > 0 ):
                   exec("year_"+typ+"_"+scenar+"_pct"+pct+"[kbasin] = year[np.min(dfp)] + int(timewin/2)")
                 else:
                   exec("year_"+typ+"_"+scenar+"_pct"+pct+"[kbasin] = 0")
                   print('No match for : '+kbasin+' '+typ+' '+scenar+' '+pct)

      if ( (basin.NAME.values[kbasin]=='LarsenC') | (basin.NAME.values[kbasin]=='Thwaites') | (basin.NAME.values[kbasin]=='LarsenB') | (basin.NAME.values[kbasin]=='Publications') | (basin.NAME.values[kbasin]=='Ronne') ):
         print('#####################################################################')
         print(basin.NAME.values[kbasin],'  SSP126 :')
         print(' ')
         print(avgru_bas_ssp126_pct05[:,kbasin])
         print(' ')
         print(np.min(avgru_bas_ssp126_pct05[:,kbasin]),np.max(avgru_bas_ssp126_pct05[:,kbasin]))
         print(year_bas_ssp126_pct05[kbasin])
         print(' ')
         print(avgru_bas_ssp126_pct95[:,kbasin])
         print(' ')
         print(np.min(avgru_bas_ssp126_pct95[:,kbasin]),np.max(avgru_bas_ssp126_pct95[:,kbasin]))
         print(year_bas_ssp126_pct95[kbasin])
         print(' ')
         print(basin.NAME.values[kbasin],'  SSP585 :')
         print(' ')
         print(avgru_bas_ssp585_pct05[:,kbasin])
         print(' ')
         print(np.min(avgru_bas_ssp585_pct05[:,kbasin]),np.max(avgru_bas_ssp585_pct05[:,kbasin]))
         print(year_bas_ssp585_pct05[kbasin])
         print(' ')
         print(avgru_bas_ssp585_pct95[:,kbasin])
         print(' ')
         print(np.min(avgru_bas_ssp585_pct95[:,kbasin]),np.max(avgru_bas_ssp585_pct95[:,kbasin]))
         print(year_bas_ssp585_pct95[kbasin])
         print(' ')

   np.savez(filez_runoff_dates,\
      threshold = threshold, \
      timewin = timewin, \
      year_isf_ssp126_pct05 = year_isf_ssp126_pct05, \
      year_isf_ssp126_pct17 = year_isf_ssp126_pct17, \
      year_isf_ssp126_pct83 = year_isf_ssp126_pct83, \
      year_isf_ssp126_pct95 = year_isf_ssp126_pct95, \
      year_isf_ssp585_pct05 = year_isf_ssp585_pct05, \
      year_isf_ssp585_pct17 = year_isf_ssp585_pct17, \
      year_isf_ssp585_pct83 = year_isf_ssp585_pct83, \
      year_isf_ssp585_pct95 = year_isf_ssp585_pct95, \
      year_bas_ssp126_pct05 = year_bas_ssp126_pct05, \
      year_bas_ssp126_pct17 = year_bas_ssp126_pct17, \
      year_bas_ssp126_pct83 = year_bas_ssp126_pct83, \
      year_bas_ssp126_pct95 = year_bas_ssp126_pct95, \
      year_bas_ssp585_pct05 = year_bas_ssp585_pct05, \
      year_bas_ssp585_pct17 = year_bas_ssp585_pct17, \
      year_bas_ssp585_pct83 = year_bas_ssp585_pct83, \
      year_bas_ssp585_pct95 = year_bas_ssp585_pct95 )

else:

   print('WARNING : starting from existing npz file : '+filez_runoff_dates)
   zz=np.load(filez_runoff_dates)
   threshold = zz['threshold']
   timewin = zz['timewin']
   print('threshold = ',threshold,'   timewin = ',timewin)
   year_isf_ssp126_pct05 = np.array(zz['year_isf_ssp126_pct05'])
   year_isf_ssp126_pct17 = np.array(zz['year_isf_ssp126_pct17'])
   year_isf_ssp126_pct83 = np.array(zz['year_isf_ssp126_pct83'])
   year_isf_ssp126_pct95 = np.array(zz['year_isf_ssp126_pct95'])
   year_isf_ssp585_pct05 = np.array(zz['year_isf_ssp585_pct05'])
   year_isf_ssp585_pct17 = np.array(zz['year_isf_ssp585_pct17'])
   year_isf_ssp585_pct83 = np.array(zz['year_isf_ssp585_pct83'])
   year_isf_ssp585_pct95 = np.array(zz['year_isf_ssp585_pct95'])
   year_bas_ssp126_pct05 = np.array(zz['year_bas_ssp126_pct05'])
   year_bas_ssp126_pct17 = np.array(zz['year_bas_ssp126_pct17'])
   year_bas_ssp126_pct83 = np.array(zz['year_bas_ssp126_pct83'])
   year_bas_ssp126_pct95 = np.array(zz['year_bas_ssp126_pct95'])
   year_bas_ssp585_pct05 = np.array(zz['year_bas_ssp585_pct05'])
   year_bas_ssp585_pct17 = np.array(zz['year_bas_ssp585_pct17'])
   year_bas_ssp585_pct83 = np.array(zz['year_bas_ssp585_pct83'])
   year_bas_ssp585_pct95 = np.array(zz['year_bas_ssp585_pct95'])

#=================================================================================

for kbasin in np.arange(Nbasin):
   if ( year_bas_ssp126_pct17[kbasin] > 0 ):
      print('==========', basin.NAME.values[kbasin], '==========')
      print('  SSP1-2.6 : ',year_bas_ssp126_pct95[kbasin], year_bas_ssp126_pct83[kbasin], year_bas_ssp126_pct17[kbasin], year_bas_ssp126_pct05[kbasin])
      print('  SSP5-8.5 : ',year_bas_ssp585_pct95[kbasin], year_bas_ssp585_pct83[kbasin], year_bas_ssp585_pct17[kbasin], year_bas_ssp585_pct05[kbasin])
