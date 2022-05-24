import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os

if not os.path.isfile('TS_profiles_ContShelf.npz'):

  print('Recalculating TS profiles (this may take a while)...')  
  
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
  
  #----------
  
  file_mshA = '/Users/njourdain/OUTPUT_AMU/mesh_mask_AMUXL12_BedMachineAntarctica-2019-05-24.nc'
  file_mshB = '/Users/njourdain/OUTPUT_AMU/mesh_mask_AMUXL12_BedMachineAntarctica-2020-07-15_v02_ICB380.nc'
  #file_mshC = '/Users/njourdain/OUTPUT_AMU/mesh_mask_AMUXL12_BedMachineAntarctica-2020-07-15_v02_ICB380.nc'
  
  ncmshA = xr.open_dataset(file_mshA,decode_cf=False)
  ncmshB = xr.open_dataset(file_mshB,decode_cf=False)
  #ncmshC = xr.open_dataset(file_mshC,decode_cf=False)
  
  file_bat = '/Users/njourdain/OUTPUT_AMU/bathy_meter_AMUXL12_BedMachineAntarctica-2020-07-15_v02.nc'
  ncbat = xr.open_dataset(file_bat,decode_cf=False)
  
  lon2d = ncmshA.glamt.isel(t=0)
  lat2d = ncmshA.gphit.isel(t=0)
  
  bathy = ncbat.Bathymetry
  
  #----------
  # T, S average over the continental shelf :
  
  mz=ncmshA.tmask.shape[1]
  print('mz = ',mz)
  mean_Tp=np.zeros((mz))*np.nan
  mean_Tf=np.zeros((mz))*np.nan
  mean_Sp=np.zeros((mz))*np.nan
  mean_Sf=np.zeros((mz))*np.nan
  
  for kk in np.arange(54):
  #for kk in np.arange(mz):
    print('   kk = ',kk)
    msk_shelfA = ncmshA.tmask.isel(t=0,z=kk).where((ncmshA.tmask.isel(t=0,z=kk)==1)&(bathy<1500)&(lon2d>=-135)&(lon2d<=-100.0)&(lat2d<-70.2),0.e0)
    msk_shelfB = ncmshB.tmask.isel(t=0,z=kk).where((ncmshB.tmask.isel(t=0,z=kk)==1)&(bathy<1500)&(lon2d>=-135)&(lon2d<=-100.0)&(lat2d<-70.2),0.e0)
    msk_shelfC = msk_shelfB
    #--
    areacella = ncmshA.e1t.isel(t=0) * ncmshA.e2t.isel(t=0)
    tmpA_area = msk_shelfA * areacella
    tmpB_area = msk_shelfB * areacella
    tmpC_area = tmpB_area
    #--
    tmp_Tp =   ncpA.toce.mean(axis=0).isel(deptht=kk) * tmpA_area \
             + ncpB.toce.mean(axis=0).isel(deptht=kk) * tmpB_area \
             + ncpC.toce.mean(axis=0).isel(deptht=kk) * tmpC_area
    tmp_Tf =   ncfA.toce.mean(axis=0).isel(deptht=kk) * tmpA_area \
             + ncfB.toce.mean(axis=0).isel(deptht=kk) * tmpB_area \
             + ncfC.toce.mean(axis=0).isel(deptht=kk) * tmpC_area
    tmp_Sp =   ncpA.soce.mean(axis=0).isel(deptht=kk) * tmpA_area \
             + ncpB.soce.mean(axis=0).isel(deptht=kk) * tmpB_area \
             + ncpC.soce.mean(axis=0).isel(deptht=kk) * tmpC_area
    tmp_Sf =   ncfA.soce.mean(axis=0).isel(deptht=kk) * tmpA_area \
             + ncfB.soce.mean(axis=0).isel(deptht=kk) * tmpB_area \
             + ncfC.soce.mean(axis=0).isel(deptht=kk) * tmpC_area
    #--
    sum_area = tmpA_area.sum(axis=(0,1)) + tmpB_area.sum(axis=(0,1)) + tmpC_area.sum(axis=(0,1))
    txp_Tp = tmp_Tp.sum(axis=(0,1)) / sum_area
    txp_Tf = tmp_Tf.sum(axis=(0,1)) / sum_area
    txp_Sp = tmp_Sp.sum(axis=(0,1)) / sum_area
    txp_Sf = tmp_Sf.sum(axis=(0,1)) / sum_area
    #--
    mean_Tp[kk] = txp_Tp.values
    mean_Tf[kk] = txp_Tf.values
    mean_Sp[kk] = txp_Sp.values
    mean_Sf[kk] = txp_Sf.values
  
  ZZ=ncmshA.gdept_1d.isel(t=0)
  
  np.savez('TS_profiles_ContShelf.npz',mean_Tp=mean_Tp,mean_Tf=mean_Tf,mean_Sp=mean_Sp,mean_Sf=mean_Sf,ZZ=ZZ)

else:

  print('  ')
  print('WARNING !!! USING EXISTING TS_profiles_ContShelf.npz !!!')
  print('  Delete before running this script if you want to regenerate it.')
  print('  ')

  npzfile = np.load('TS_profiles_ContShelf.npz')
  mean_Tp = npzfile['mean_Tp']
  mean_Tf = npzfile['mean_Tf']
  mean_Sp = npzfile['mean_Sp']
  mean_Sf = npzfile['mean_Sf']
  ZZ = npzfile['ZZ']

#----------

fig, axs = plt.subplots(nrows=1,ncols=2,figsize=(21.0,7.0))
axs = axs.ravel()

colp='cornflowerblue'
colf='orange'

lineTp, = axs[0].plot(mean_Tp,ZZ,linewidth=1.5,color=colp)
lineTf, = axs[0].plot(mean_Tf,ZZ,linewidth=1.5,color=colf)
axs[0].set_xlabel('(a) Conservative Temperature (degC)',fontsize=16)
axs[0].set_ylabel('Depth (m)',fontsize=14)
axs[0].xaxis.set_ticks_position('top')
axs[0].xaxis.set_label_position('top')
axs[0].invert_yaxis()
axs[0].set_ylim(1250,0)
axs[0].tick_params(axis='x', labelsize=12)
axs[0].tick_params(axis='y', labelsize=12)
axs[0].grid(True) 
axs[0].legend((lineTp,lineTf),('1989-2009, mean(A,B,C)','2080-2100, mean(A,B,C)'),fontsize=14)

lineSp, = axs[1].plot(mean_Sp,ZZ,linewidth=1.5,color=colp)
lineSf, = axs[1].plot(mean_Sf,ZZ,linewidth=1.5,color=colf)
axs[1].set_xlabel('(b) Absolute Salinity (g/kg)',fontsize=16)
axs[1].set_ylabel('Depth (m)',fontsize=14)
axs[1].xaxis.set_ticks_position('top')
axs[1].xaxis.set_label_position('top')
axs[1].invert_yaxis()
axs[1].set_ylim(1250,0)
axs[1].tick_params(axis='x', labelsize=12)
axs[1].tick_params(axis='y', labelsize=12)
axs[1].grid(True) 
axs[1].legend((lineSp,lineSf),('1989-2009, mean(A,B,C)','2080-2100, mean(A,B,C)'),fontsize=14)

fig.savefig('TS_profiles_ContShelf.jpg')
fig.savefig('TS_profiles_ContShelf.pdf')
