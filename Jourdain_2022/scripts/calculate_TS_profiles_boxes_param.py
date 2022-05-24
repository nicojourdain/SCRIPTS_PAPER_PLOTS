import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os
import functions_nico as nico

RT=6370. # Earth radius [km]

file_mshA = '/Users/njourdain/OUTPUT_AMU/mesh_mask_AMUXL12_BedMachineAntarctica-2019-05-24.nc'
file_mshB = '/Users/njourdain/OUTPUT_AMU/mesh_mask_AMUXL12_BedMachineAntarctica-2020-07-15_v02_ICB380.nc'
#file_mshC = '/Users/njourdain/OUTPUT_AMU/mesh_mask_AMUXL12_BedMachineAntarctica-2020-07-15_v02_ICB380.nc'

ncmshA = xr.open_dataset(file_mshA,decode_cf=False)
ncmshB = xr.open_dataset(file_mshB,decode_cf=False)
#ncmshC = xr.open_dataset(file_mshC,decode_cf=False)

file_isfA = '/Users/njourdain/OUTPUT_AMU/isf_mask_AMUXL12_BedMachineAntarctica-2019-05-24.nc'
file_isfB = '/Users/njourdain/OUTPUT_AMU/isf_mask_AMUXL12_BedMachineAntarctica-2020-07-15_v02.nc'
#file_isfC = '/Users/njourdain/OUTPUT_AMU/isf_mask_AMUXL12_BedMachineAntarctica-2020-07-15_v02.nc'

ncisfA = xr.open_dataset(file_isfA,decode_cf=False,concat_characters=True)
ncisfB = xr.open_dataset(file_isfB,decode_cf=False,concat_characters=True)
#ncisfC = xr.open_dataset(file_isfC,decode_cf=False)

#isflist = np.array(['Getz', 'Dotson', 'Crosson', 'Thwaites', 'Pine Island', 'Cosgrove', 'Abbot', 'Venable' ], dtype=str)
#isfnum  = np.array([ 22-1 ,  52-1   ,  38-1    ,  23-1     ,  66-1        ,  53-1     ,  39-1  ,  24-1     ], dtype=int)
isflist = np.array(['Getz', 'Dotson', 'Crosson', 'Thwaites', 'Pine Island', 'Cosgrove' ], dtype=str)
isfnum  = np.array([ 22-1 ,  52-1   ,  38-1    ,  23-1     ,  66-1        ,  53-1      ], dtype=int)
#isflist = np.array([ 'Pine Island' ], dtype=str)
#isfnum  = np.array([  66-1         ], dtype=int)

Ncav = np.size(isfnum)
mz = ncmshA.tmask.shape[1]
areacella = ncmshA.e1t.isel(t=0) * ncmshA.e2t.isel(t=0)
ZZ=ncmshA.gdept_1d.isel(t=0)

#====================================================================
# T, S data

file_TS_p_A = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MAR/climato_monthly_AMUXL12-GNJ002_BM02MAR_grid_T_1989_2009.nc'
file_TS_f_A = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MARrcp85/climato_monthly_AMUXL12-GNJ002_BM02MARrcp85_grid_T_1989_2009.nc'
file_TS_p_B = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM03MAR/climato_monthly_AMUXL12-GNJ002_BM03MAR_grid_T_1989_2009.nc'
file_TS_f_B = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM03MARrcBDY/climato_monthly_AMUXL12-GNJ002_BM03MARrcBDY_grid_T_1989_2009.nc'
file_TS_p_C = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM04MAR/climato_monthly_AMUXL12-GNJ002_BM04MAR_grid_T_1989_2009.nc'
file_TS_f_C = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM04MARrcp85/climato_monthly_AMUXL12-GNJ002_BM04MARrcp85_grid_T_1989_2009.nc'
  
ncpA = xr.open_dataset(file_TS_p_A,decode_cf=False)
ncpB = xr.open_dataset(file_TS_p_B,decode_cf=False)
ncpC = xr.open_dataset(file_TS_p_C,decode_cf=False)
ncfA = xr.open_dataset(file_TS_f_A,decode_cf=False)
ncfB = xr.open_dataset(file_TS_f_B,decode_cf=False)
ncfC = xr.open_dataset(file_TS_f_C,decode_cf=False)

#====================================================================
# Calculate T,S profiles for all ice shelves:
# (using same boxes for A, B, C)

mean_TpA = np.zeros((Ncav,mz)) * np.nan
mean_TfA = np.zeros((Ncav,mz)) * np.nan
mean_SpA = np.zeros((Ncav,mz)) * np.nan
mean_SfA = np.zeros((Ncav,mz)) * np.nan
mean_TpB = np.zeros((Ncav,mz)) * np.nan
mean_TfB = np.zeros((Ncav,mz)) * np.nan
mean_SpB = np.zeros((Ncav,mz)) * np.nan
mean_SfB = np.zeros((Ncav,mz)) * np.nan
mean_TpC = np.zeros((Ncav,mz)) * np.nan
mean_TfC = np.zeros((Ncav,mz)) * np.nan
mean_SpC = np.zeros((Ncav,mz)) * np.nan
mean_SfC = np.zeros((Ncav,mz)) * np.nan

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

   mskbox = ncmshB.tmask.isel(t=0,z=0).where(   (ncmshB.glamt.isel(t=0)<lonmax+dlon) \
                                              & (ncmshB.glamt.isel(t=0)>lonmin-dlon) \
                                              & (ncmshB.gphit.isel(t=0)<latmax+dlat) \
                                              & (ncmshB.gphit.isel(t=0)>latmin-dlat) \
                                             , 0.e0 )

   for kk in np.arange(60):
   #for kk in np.arange(mz):
 
     print('   kk = ',kk)
     #--
     tmpA_area = ncmshA.tmask.isel(t=0,z=kk) * mskbox * areacella
     tmpB_area = ncmshB.tmask.isel(t=0,z=kk) * mskbox * areacella
     tmpC_area = tmpB_area
     #-----------
     tmp_TpA =  ncpA.toce.mean(axis=0).isel(deptht=kk) * tmpA_area
     tmp_TpB =  ncpB.toce.mean(axis=0).isel(deptht=kk) * tmpB_area
     tmp_TpC =  ncpC.toce.mean(axis=0).isel(deptht=kk) * tmpC_area
     #--
     tmp_TfA =  ncfA.toce.mean(axis=0).isel(deptht=kk) * tmpA_area
     tmp_TfB =  ncfB.toce.mean(axis=0).isel(deptht=kk) * tmpB_area
     tmp_TfC =  ncfC.toce.mean(axis=0).isel(deptht=kk) * tmpC_area
     #--
     tmp_SpA =  ncpA.soce.mean(axis=0).isel(deptht=kk) * tmpA_area
     tmp_SpB =  ncpB.soce.mean(axis=0).isel(deptht=kk) * tmpB_area
     tmp_SpC =  ncpC.soce.mean(axis=0).isel(deptht=kk) * tmpC_area
     #--
     tmp_SfA =  ncfA.soce.mean(axis=0).isel(deptht=kk) * tmpA_area
     tmp_SfB =  ncfB.soce.mean(axis=0).isel(deptht=kk) * tmpB_area
     tmp_SfC =  ncfC.soce.mean(axis=0).isel(deptht=kk) * tmpC_area
     #--
     sum_areaA = tmpA_area.sum(axis=(0,1))
     sum_areaB = tmpB_area.sum(axis=(0,1))
     sum_areaC = sum_areaB
     #----------
     txp_TpA = tmp_TpA.sum(axis=(0,1)) / sum_areaA
     txp_TpB = tmp_TpB.sum(axis=(0,1)) / sum_areaB
     txp_TpC = tmp_TpC.sum(axis=(0,1)) / sum_areaC
     #--
     txp_TfA = tmp_TfA.sum(axis=(0,1)) / sum_areaA
     txp_TfB = tmp_TfB.sum(axis=(0,1)) / sum_areaB
     txp_TfC = tmp_TfC.sum(axis=(0,1)) / sum_areaC
     #--
     txp_SpA = tmp_SpA.sum(axis=(0,1)) / sum_areaA
     txp_SpB = tmp_SpB.sum(axis=(0,1)) / sum_areaB
     txp_SpC = tmp_SpC.sum(axis=(0,1)) / sum_areaC
     #--
     txp_SfA = tmp_SfA.sum(axis=(0,1)) / sum_areaA
     txp_SfB = tmp_SfB.sum(axis=(0,1)) / sum_areaB
     txp_SfC = tmp_SfC.sum(axis=(0,1)) / sum_areaC
     #----------
     mean_TpA[kcav,kk] = txp_TpA.values
     mean_TpB[kcav,kk] = txp_TpB.values
     mean_TpC[kcav,kk] = txp_TpC.values
     #--
     mean_TfA[kcav,kk] = txp_TfA.values
     mean_TfB[kcav,kk] = txp_TfB.values
     mean_TfC[kcav,kk] = txp_TfC.values
     #--
     mean_SpA[kcav,kk] = txp_SpA.values
     mean_SpB[kcav,kk] = txp_SpB.values
     mean_SpC[kcav,kk] = txp_SpC.values
     #--
     mean_SfA[kcav,kk] = txp_SfA.values
     mean_SfB[kcav,kk] = txp_SfB.values
     mean_SfC[kcav,kk] = txp_SfC.values
     #--
 
np.savez('TS_profiles_boxes_param.npz',mean_TpA=mean_TpA,\
                                       mean_TpB=mean_TpB,\
                                       mean_TpC=mean_TpC,\
                                       mean_TfA=mean_TfA,\
                                       mean_TfB=mean_TfB,\
                                       mean_TfC=mean_TfC,\
                                       mean_SpA=mean_SpA,\
                                       mean_SpB=mean_SpB,\
                                       mean_SpC=mean_SpC,\
                                       mean_SfA=mean_SfA,\
                                       mean_SfB=mean_SfB,\
                                       mean_SfC=mean_SfC,\
                                       isflist=isflist,\
                                       isfnum=isfnum,\
                                       ZZ=ZZ)
