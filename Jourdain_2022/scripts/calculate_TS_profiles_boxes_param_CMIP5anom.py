import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os
import functions_nico as nico

RT=6370. # Earth radius [km]

file_zmsh = 'CMIP5_DATA/ORCA025.L75-MJM91_mesh_zgr_REDUCED.nc'
nczmsh = xr.open_dataset(file_zmsh,decode_cf=False)

file_hmsh = 'CMIP5_DATA/ORCA025.L75-MJM91_mesh_hgr_REDUCED.nc'
nchmsh = xr.open_dataset(file_hmsh,decode_cf=False)

file_mask = 'CMIP5_DATA/ORCA025.L75-MJM91_byte_mask_REDUCED.nc'
ncmask = xr.open_dataset(file_mask,decode_cf=False)

file_isfA = '/Users/njourdain/OUTPUT_AMU/isf_mask_AMUXL12_BedMachineAntarctica-2019-05-24.nc'
ncisfA = xr.open_dataset(file_isfA,decode_cf=False,concat_characters=True)

#isflist = np.array(['Getz', 'Dotson', 'Crosson', 'Thwaites', 'Pine Island', 'Cosgrove', 'Abbot', 'Venable' ], dtype=str)
#isfnum  = np.array([ 22-1 ,  52-1   ,  38-1    ,  23-1     ,  66-1        ,  53-1     ,  39-1  ,  24-1     ], dtype=int)
isflist = np.array(['Getz', 'Dotson', 'Crosson', 'Thwaites', 'Pine Island', 'Cosgrove' ], dtype=str)
isfnum  = np.array([ 22-1 ,  52-1   ,  38-1    ,  23-1     ,  66-1        ,  53-1      ], dtype=int)
#isflist = np.array([ 'Pine Island' ], dtype=str)
#isfnum  = np.array([  66-1         ], dtype=int)

Ncav = np.size(isfnum)
mz = ncmask.tmask.shape[1]
areacella = nchmsh.e1t.isel(t=0) * nchmsh.e2t.isel(t=0)
ZZ=nczmsh.gdept_0.isel(t=0)

#====================================================================
# T, S data

file_T = 'CMIP5_DATA/Oce_anom_ALL31_CMIP5_votemper_gridT_REDUCED.nc'
file_S = 'CMIP5_DATA/Oce_anom_ALL31_CMIP5_vosaline_gridT_REDUCED.nc'
  
ncT = xr.open_dataset(file_T,decode_cf=False)
ncS = xr.open_dataset(file_S,decode_cf=False)

#====================================================================
# Calculate T,S profiles for all ice shelves:
# (using same boxes for A, B, C)

mean_Tanom = np.zeros((Ncav,mz)) * np.nan
mean_Sanom = np.zeros((Ncav,mz)) * np.nan

for kcav in np.arange(Ncav):

   print(isfnum[kcav], isflist[kcav])

   lonmax = ncisfA.front_max_lon.isel(Nisf=isfnum[kcav]) 
   lonmin = ncisfA.front_min_lon.isel(Nisf=isfnum[kcav])
   latmax = ncisfA.front_max_lat.isel(Nisf=isfnum[kcav])  
   latmin = ncisfA.front_min_lat.isel(Nisf=isfnum[kcav])
   print(lonmax.values,lonmin.values,latmax.values,latmin.values)

   # delta lon,lat corresponding to ~50km :
   dlon = 50. * 180. / ( np.pi * RT * np.cos(0.5*(latmin+latmax)*np.pi/180.) )
   dlat = 50. * 180. / ( np.pi * RT )

   mskbox = ncmask.tmask.isel(t=0,z=0).where(   (nchmsh.glamt.isel(t=0)<lonmax+dlon) \
                                              & (nchmsh.glamt.isel(t=0)>lonmin-dlon) \
                                              & (nchmsh.gphit.isel(t=0)<latmax+dlat) \
                                              & (nchmsh.gphit.isel(t=0)>latmin-dlat) \
                                             , 0.e0 )

   for kk in np.arange(60):
   #for kk in np.arange(mz):
 
     print('   kk = ',kk)
     #--
     tmp_area = ncmask.tmask.isel(t=0,z=kk) * mskbox * areacella
     #-----------
     tmp_T =  ncT.votemper_anom.mean(axis=0).isel(z=kk) * tmp_area
     tmp_S =  ncS.vosaline_anom.mean(axis=0).isel(z=kk) * tmp_area
     #--
     sum_area = tmp_area.sum(axis=(0,1))
     #----------
     txp_T = tmp_T.sum(axis=(0,1)) / sum_area
     txp_S = tmp_S.sum(axis=(0,1)) / sum_area
     #----------
     mean_Tanom[kcav,kk] = txp_T.values
     mean_Sanom[kcav,kk] = txp_S.values
     #--
 
np.savez('TS_profiles_boxes_param_CMIP5anom.npz',mean_Tanom=mean_Tanom,\
                                                 mean_Sanom=mean_Sanom,\
                                                 isflist=isflist,\
                                                 isfnum=isfnum,\
                                                 ZZ=ZZ)

print(mean_Tanom)
