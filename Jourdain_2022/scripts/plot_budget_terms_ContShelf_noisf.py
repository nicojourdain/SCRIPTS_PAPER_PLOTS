import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os
import functions_nico as nico

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

#=========================================================================================
# dT terms
#=========================================================================================

npz_TS_file = np.load('TS_profiles_ContShelf_noisf.npz')
mean_dT = npz_TS_file['mean_Tf'] - npz_TS_file['mean_Tp']

if not os.path.isfile('budget_dT_terms_ContShelf_noisf.npz'):

  print('Recalculating profiles of heat and salt budget terms (this may take a while)...')  
  
  file_T_p_A = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MAR/deltaT_ttrd_climato_AMUXL12-GNJ002_BM02MAR_since_1979_climato_1989-2009.nc'
  file_T_f_A = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MARrcp85/deltaT_ttrd_climato_AMUXL12-GNJ002_BM02MARrcp85_since_1979_climato_1989-2009.nc'
  file_T_p_B = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM03MAR/deltaT_ttrd_climato_AMUXL12-GNJ002_BM03MAR_since_1979_climato_1989-2009.nc'
  file_T_f_B = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM03MARrcBDY/deltaT_ttrd_climato_AMUXL12-GNJ002_BM03MARrcBDY_since_1979_climato_1989-2009.nc'
  file_T_p_C = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM04MAR/deltaT_ttrd_climato_AMUXL12-GNJ002_BM04MAR_since_1979_climato_1989-2009.nc'
  file_T_f_C = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM04MARrcp85/deltaT_ttrd_climato_AMUXL12-GNJ002_BM04MARrcp85_since_1979_climato_1989-2009.nc'

  ncTpA = xr.open_dataset(file_T_p_A,decode_cf=False)
  ncTfA = xr.open_dataset(file_T_f_A,decode_cf=False)
  ncTpB = xr.open_dataset(file_T_p_B,decode_cf=False)
  ncTfB = xr.open_dataset(file_T_f_B,decode_cf=False)
  ncTpC = xr.open_dataset(file_T_p_C,decode_cf=False)
  ncTfC = xr.open_dataset(file_T_f_C,decode_cf=False)

  #----------
  # Average over the continental shelf :
  
  mz=ncmshA.tmask.shape[1]
  print('mz = ',mz)
  mean_dT_hadv=np.zeros((mz))*np.nan
  mean_dT_vadv=np.zeros((mz))*np.nan
  mean_dT_hmix=np.zeros((mz))*np.nan
  mean_dT_qsol=np.zeros((mz))*np.nan
  mean_dT_vmxq=np.zeros((mz))*np.nan
  mean_dT_qisf=np.zeros((mz))*np.nan
  
  for kk in np.arange(54):
  #for kk in np.arange(mz):
    print('   kk = ',kk)
    msk_shelfA = ncmshA.tmask.isel(t=0,z=kk).where((ncmshA.tmask.isel(t=0,z=kk)==1)&(ncmshA.misf.isel(t=0)==1)&(bathy<1500)&(lon2d>=-135)&(lon2d<=-100.0)&(lat2d<-70.2),0.e0)
    msk_shelfB = ncmshB.tmask.isel(t=0,z=kk).where((ncmshB.tmask.isel(t=0,z=kk)==1)&(ncmshB.misf.isel(t=0)==1)&(bathy<1500)&(lon2d>=-135)&(lon2d<=-100.0)&(lat2d<-70.2),0.e0)
    msk_shelfC = msk_shelfB
    #--
    areacella = ncmshA.e1t.isel(t=0) * ncmshA.e2t.isel(t=0)
    tmpA_area = msk_shelfA * areacella
    tmpB_area = msk_shelfB * areacella
    tmpC_area = tmpB_area
    #--
    # Horizontal advection:
    tmp_dT_hadv =   ( ncTfA.dT___xad.isel(deptht=kk) - ncTpA.dT___xad.isel(deptht=kk) ) * tmpA_area \
                  + ( ncTfB.dT___xad.isel(deptht=kk) - ncTpB.dT___xad.isel(deptht=kk) ) * tmpB_area \
                  + ( ncTfC.dT___xad.isel(deptht=kk) - ncTpC.dT___xad.isel(deptht=kk) ) * tmpC_area \
                  + ( ncTfA.dT___yad.isel(deptht=kk) - ncTpA.dT___yad.isel(deptht=kk) ) * tmpA_area \
                  + ( ncTfB.dT___yad.isel(deptht=kk) - ncTpB.dT___yad.isel(deptht=kk) ) * tmpB_area \
                  + ( ncTfC.dT___yad.isel(deptht=kk) - ncTpC.dT___yad.isel(deptht=kk) ) * tmpC_area
    # Vertical advection:
    tmp_dT_vadv =   ( ncTfA.dT___zad.isel(deptht=kk) - ncTpA.dT___zad.isel(deptht=kk) ) * tmpA_area \
                  + ( ncTfB.dT___zad.isel(deptht=kk) - ncTpB.dT___zad.isel(deptht=kk) ) * tmpB_area \
                  + ( ncTfC.dT___zad.isel(deptht=kk) - ncTpC.dT___zad.isel(deptht=kk) ) * tmpC_area
    # Horizontal mixing:
    tmp_dT_hmix =   ( ncTfA.dT___ldf.isel(deptht=kk) - ncTpA.dT___ldf.isel(deptht=kk) ) * tmpA_area \
                  + ( ncTfB.dT___ldf.isel(deptht=kk) - ncTpB.dT___ldf.isel(deptht=kk) ) * tmpB_area \
                  + ( ncTfC.dT___ldf.isel(deptht=kk) - ncTpC.dT___ldf.isel(deptht=kk) ) * tmpC_area
    # Solar heat flux:
    tmp_dT_qsol =   ( ncTfA.dT___qsr.isel(deptht=kk) - ncTpA.dT___qsr.isel(deptht=kk) ) * tmpA_area \
                  + ( ncTfB.dT___qsr.isel(deptht=kk) - ncTpB.dT___qsr.isel(deptht=kk) ) * tmpB_area \
                  + ( ncTfC.dT___qsr.isel(deptht=kk) - ncTpC.dT___qsr.isel(deptht=kk) ) * tmpC_area
    # Vertical mixing + surface fluxes (all except ice-shelf):
    tmp_dT_vmxq =   ( ncTfA.dT___zdf.isel(deptht=kk) - ncTpA.dT___zdf.isel(deptht=kk) ) * tmpA_area \
                  + ( ncTfB.dT___zdf.isel(deptht=kk) - ncTpB.dT___zdf.isel(deptht=kk) ) * tmpB_area \
                  + ( ncTfC.dT___zdf.isel(deptht=kk) - ncTpC.dT___zdf.isel(deptht=kk) ) * tmpC_area
    if kk == 0:
      tmp_dT_vmxq =   tmp_dT_vmxq \
                    + ( ncTfA.dT___qns.isel(deptht=kk) - ncTpA.dT___qns.isel(deptht=kk) ) * tmpA_area \
                    + ( ncTfB.dT___qns.isel(deptht=kk) - ncTpB.dT___qns.isel(deptht=kk) ) * tmpB_area \
                    + ( ncTfC.dT___qns.isel(deptht=kk) - ncTpC.dT___qns.isel(deptht=kk) ) * tmpC_area
    # Ice-shelf (all non-solar at depth):
    tmp_dT_qisf =   ( ncTfA.dT___qns.isel(deptht=kk) - ncTpA.dT___qns.isel(deptht=kk) ) * tmpA_area \
                  + ( ncTfB.dT___qns.isel(deptht=kk) - ncTpB.dT___qns.isel(deptht=kk) ) * tmpB_area \
                  + ( ncTfC.dT___qns.isel(deptht=kk) - ncTpC.dT___qns.isel(deptht=kk) ) * tmpC_area
    if kk == 0:
       tmp_dT_qisf = tmp_dT_qisf * 0.e0
    #--
    sum_area = tmpA_area.sum(axis=(0,1)) + tmpB_area.sum(axis=(0,1)) + tmpC_area.sum(axis=(0,1))
    txp_dT_hadv = tmp_dT_hadv.sum(axis=(0,1)) / sum_area
    txp_dT_vadv = tmp_dT_vadv.sum(axis=(0,1)) / sum_area
    txp_dT_hmix = tmp_dT_hmix.sum(axis=(0,1)) / sum_area
    txp_dT_qsol = tmp_dT_qsol.sum(axis=(0,1)) / sum_area
    txp_dT_vmxq = tmp_dT_vmxq.sum(axis=(0,1)) / sum_area
    txp_dT_qisf = tmp_dT_qisf.sum(axis=(0,1)) / sum_area
    #--
    mean_dT_hadv[kk] = txp_dT_hadv.values
    mean_dT_vadv[kk] = txp_dT_vadv.values
    mean_dT_hmix[kk] = txp_dT_hmix.values
    mean_dT_qsol[kk] = txp_dT_qsol.values
    mean_dT_vmxq[kk] = txp_dT_vmxq.values
    mean_dT_qisf[kk] = txp_dT_qisf.values
  
  ZZ=ncmshA.gdept_1d.isel(t=0)
  
  np.savez('budget_dT_terms_ContShelf_noisf.npz',mean_dT_hadv=mean_dT_hadv,\
                                          mean_dT_vadv=mean_dT_vadv,\
                                          mean_dT_hmix=mean_dT_hmix,\
                                          mean_dT_qsol=mean_dT_qsol,\
                                          mean_dT_vmxq=mean_dT_vmxq,\
                                          mean_dT_qisf=mean_dT_qisf,\
                                          ZZ=ZZ)

else:

  print('  ')
  print('WARNING !!! USING EXISTING budget_dT_terms_ContShelf_noisf.npz !!!')
  print('  Delete before running this script if you want to regenerate it.')
  print('  ')

  npzfileT = np.load('budget_dT_terms_ContShelf_noisf.npz')
  mean_dT_hadv = npzfileT['mean_dT_hadv']
  mean_dT_vadv = npzfileT['mean_dT_vadv']
  mean_dT_hmix = npzfileT['mean_dT_hmix']
  mean_dT_qsol = npzfileT['mean_dT_qsol']
  mean_dT_vmxq = npzfileT['mean_dT_vmxq']
  mean_dT_qisf = npzfileT['mean_dT_qisf']
  ZZ = npzfileT['ZZ']

#=========================================================================================
# dS terms
#=========================================================================================

mean_dS = npz_TS_file['mean_Sf'] - npz_TS_file['mean_Sp']

if not os.path.isfile('budget_dS_terms_ContShelf_noisf.npz'):

  print('Recalculating profiles of heat and salt budget terms (this may take a while)...')  
  
  file_S_p_A = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MAR/deltaS_strd_climato_AMUXL12-GNJ002_BM02MAR_since_1979_climato_1989-2009.nc'
  file_S_f_A = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MARrcp85/deltaS_strd_climato_AMUXL12-GNJ002_BM02MARrcp85_since_1979_climato_1989-2009.nc'
  file_S_p_B = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM03MAR/deltaS_strd_climato_AMUXL12-GNJ002_BM03MAR_since_1979_climato_1989-2009.nc'
  file_S_f_B = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM03MARrcBDY/deltaS_strd_climato_AMUXL12-GNJ002_BM03MARrcBDY_since_1979_climato_1989-2009.nc'
  file_S_p_C = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM04MAR/deltaS_strd_climato_AMUXL12-GNJ002_BM04MAR_since_1979_climato_1989-2009.nc'
  file_S_f_C = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM04MARrcp85/deltaS_strd_climato_AMUXL12-GNJ002_BM04MARrcp85_since_1979_climato_1989-2009.nc'

  ncSpA = xr.open_dataset(file_S_p_A,decode_cf=False)
  ncSfA = xr.open_dataset(file_S_f_A,decode_cf=False)
  ncSpB = xr.open_dataset(file_S_p_B,decode_cf=False)
  ncSfB = xr.open_dataset(file_S_f_B,decode_cf=False)
  ncSpC = xr.open_dataset(file_S_p_C,decode_cf=False)
  ncSfC = xr.open_dataset(file_S_f_C,decode_cf=False)

  #----------
  # Average over the continental shelf :
  
  mz=ncmshA.tmask.shape[1]
  print('mz = ',mz)
  mean_dS_hadv=np.zeros((mz))*np.nan
  mean_dS_vadv=np.zeros((mz))*np.nan
  mean_dS_hmix=np.zeros((mz))*np.nan
  mean_dS_vmxq=np.zeros((mz))*np.nan
  mean_dS_qisf=np.zeros((mz))*np.nan
  
  for kk in np.arange(54):
  #for kk in np.arange(mz):
    print('   kk = ',kk)
    msk_shelfA = ncmshA.tmask.isel(t=0,z=kk).where((ncmshA.tmask.isel(t=0,z=kk)==1)&(ncmshA.misf.isel(t=0)==1)&(bathy<1500)&(lon2d>=-135)&(lon2d<=-100.0)&(lat2d<-70.2),0.e0)
    msk_shelfB = ncmshB.tmask.isel(t=0,z=kk).where((ncmshB.tmask.isel(t=0,z=kk)==1)&(ncmshB.misf.isel(t=0)==1)&(bathy<1500)&(lon2d>=-135)&(lon2d<=-100.0)&(lat2d<-70.2),0.e0)
    msk_shelfC = msk_shelfB
    #--
    areacella = ncmshA.e1t.isel(t=0) * ncmshA.e2t.isel(t=0)
    tmpA_area = msk_shelfA * areacella
    tmpB_area = msk_shelfB * areacella
    tmpC_area = tmpB_area
    #--
    # Horizontal advection:
    tmp_dS_hadv =   ( ncSfA.dS___xad.isel(deptht=kk) - ncSpA.dS___xad.isel(deptht=kk) ) * tmpA_area \
                  + ( ncSfB.dS___xad.isel(deptht=kk) - ncSpB.dS___xad.isel(deptht=kk) ) * tmpB_area \
                  + ( ncSfC.dS___xad.isel(deptht=kk) - ncSpC.dS___xad.isel(deptht=kk) ) * tmpC_area \
                  + ( ncSfA.dS___yad.isel(deptht=kk) - ncSpA.dS___yad.isel(deptht=kk) ) * tmpA_area \
                  + ( ncSfB.dS___yad.isel(deptht=kk) - ncSpB.dS___yad.isel(deptht=kk) ) * tmpB_area \
                  + ( ncSfC.dS___yad.isel(deptht=kk) - ncSpC.dS___yad.isel(deptht=kk) ) * tmpC_area
    # Vertical advection:
    tmp_dS_vadv =   ( ncSfA.dS___zad.isel(deptht=kk) - ncSpA.dS___zad.isel(deptht=kk) ) * tmpA_area \
                  + ( ncSfB.dS___zad.isel(deptht=kk) - ncSpB.dS___zad.isel(deptht=kk) ) * tmpB_area \
                  + ( ncSfC.dS___zad.isel(deptht=kk) - ncSpC.dS___zad.isel(deptht=kk) ) * tmpC_area
    # Horizontal mixing:
    tmp_dS_hmix =   ( ncSfA.dS___ldf.isel(deptht=kk) - ncSpA.dS___ldf.isel(deptht=kk) ) * tmpA_area \
                  + ( ncSfB.dS___ldf.isel(deptht=kk) - ncSpB.dS___ldf.isel(deptht=kk) ) * tmpB_area \
                  + ( ncSfC.dS___ldf.isel(deptht=kk) - ncSpC.dS___ldf.isel(deptht=kk) ) * tmpC_area
    # Vertical mixing + surface fluxes (all except ice-shelf):
    tmp_dS_vmxq =   ( ncSfA.dS___zdf.isel(deptht=kk) - ncSpA.dS___zdf.isel(deptht=kk) ) * tmpA_area \
                  + ( ncSfB.dS___zdf.isel(deptht=kk) - ncSpB.dS___zdf.isel(deptht=kk) ) * tmpB_area \
                  + ( ncSfC.dS___zdf.isel(deptht=kk) - ncSpC.dS___zdf.isel(deptht=kk) ) * tmpC_area
    if kk == 0:
      tmp_dS_vmxq =   tmp_dS_vmxq \
                    + ( ncSfA.dS___cdt.isel(deptht=kk) - ncSpA.dS___cdt.isel(deptht=kk) ) * tmpA_area \
                    + ( ncSfB.dS___cdt.isel(deptht=kk) - ncSpB.dS___cdt.isel(deptht=kk) ) * tmpB_area \
                    + ( ncSfC.dS___cdt.isel(deptht=kk) - ncSpC.dS___cdt.isel(deptht=kk) ) * tmpC_area
    # Ice-shelf (all non-solar at depth):
    tmp_dS_qisf =   ( ncSfA.dS___cdt.isel(deptht=kk) - ncSpA.dS___cdt.isel(deptht=kk) ) * tmpA_area \
                  + ( ncSfB.dS___cdt.isel(deptht=kk) - ncSpB.dS___cdt.isel(deptht=kk) ) * tmpB_area \
                  + ( ncSfC.dS___cdt.isel(deptht=kk) - ncSpC.dS___cdt.isel(deptht=kk) ) * tmpC_area
    if kk == 0:
       tmp_dS_qisf = tmp_dS_qisf * 0.e0
    #--
    sum_area = tmpA_area.sum(axis=(0,1)) + tmpB_area.sum(axis=(0,1)) + tmpC_area.sum(axis=(0,1))
    txp_dS_hadv = tmp_dS_hadv.sum(axis=(0,1)) / sum_area
    txp_dS_vadv = tmp_dS_vadv.sum(axis=(0,1)) / sum_area
    txp_dS_hmix = tmp_dS_hmix.sum(axis=(0,1)) / sum_area
    txp_dS_vmxq = tmp_dS_vmxq.sum(axis=(0,1)) / sum_area
    txp_dS_qisf = tmp_dS_qisf.sum(axis=(0,1)) / sum_area
    #--
    mean_dS_hadv[kk] = txp_dS_hadv.values
    mean_dS_vadv[kk] = txp_dS_vadv.values
    mean_dS_hmix[kk] = txp_dS_hmix.values
    mean_dS_vmxq[kk] = txp_dS_vmxq.values
    mean_dS_qisf[kk] = txp_dS_qisf.values
  
  ZZ=ncmshA.gdept_1d.isel(t=0)
  
  np.savez('budget_dS_terms_ContShelf_noisf.npz',mean_dS_hadv=mean_dS_hadv,\
                                          mean_dS_vadv=mean_dS_vadv,\
                                          mean_dS_hmix=mean_dS_hmix,\
                                          mean_dS_vmxq=mean_dS_vmxq,\
                                          mean_dS_qisf=mean_dS_qisf,\
                                          ZZ=ZZ)

else:

  print('  ')
  print('WARNING !!! USING EXISTING budget_dS_terms_ContShelf_noisf.npz !!!')
  print('  Delete before running this script if you want to regenerate it.')
  print('  ')

  npzfileS = np.load('budget_dS_terms_ContShelf_noisf.npz')
  mean_dS_hadv = npzfileS['mean_dS_hadv']
  mean_dS_vadv = npzfileS['mean_dS_vadv']
  mean_dS_hmix = npzfileS['mean_dS_hmix']
  mean_dS_vmxq = npzfileS['mean_dS_vmxq']
  mean_dS_qisf = npzfileS['mean_dS_qisf']
  ZZ = npzfileS['ZZ']

#==========================================================================
# bug fix (cdt for isf contribution not saved):

mean_dS_qisf = mean_dS - mean_dS_hadv - mean_dS_vadv - mean_dS_hmix - mean_dS_vmxq

#==========================================================================
# PLOT :

fig, axs = plt.subplots(nrows=1,ncols=2,figsize=(21.0,7.0))
axs = axs.ravel()

col=nico.cbcol(8,withblack=True)

line11, = axs[0].plot(mean_dT,ZZ,color=col[0])
line12, = axs[0].plot(mean_dT_hadv,ZZ,color=col[1])
line13, = axs[0].plot(mean_dT_vadv,ZZ,color=col[2])
line14, = axs[0].plot(mean_dT_hmix,ZZ,color=col[3])
line15, = axs[0].plot(mean_dT_vmxq+mean_dT_qsol,ZZ,color=col[4])
line16, = axs[0].plot(mean_dT_qisf,ZZ,color=col[7])
axs[0].set_xlabel('(e) Conservative Temperature anomaly (degC)',fontsize=16)
axs[0].set_ylabel('Depth (m)',fontsize=14)
axs[0].xaxis.set_ticks_position('top')
axs[0].xaxis.set_label_position('top')
axs[0].invert_yaxis()
axs[0].set_xlim(-14,25)
axs[0].set_ylim(1250,0)
axs[0].tick_params(axis='x', labelsize=12)
axs[0].tick_params(axis='y', labelsize=12)
axs[0].grid(True) 
axs[0].text(6,700,'Without ice-shelf cavities',color='k',fontsize=14,fontweight='bold',HorizontalAlignment='left')
axs[0].legend((line11,line12,line13,line14,line15,line16),\
              ('$\Delta$ T','Horizontal Advection','Vertical Advection','Horizontal mixing','Vert. Mix. + Surf. Fluxes','Ice Shelf heat flux'),fontsize=14)

line21, = axs[1].plot(mean_dS,ZZ,color=col[0])
line22, = axs[1].plot(mean_dS_hadv,ZZ,color=col[1])
line23, = axs[1].plot(mean_dS_vadv,ZZ,color=col[2])
line24, = axs[1].plot(mean_dS_hmix,ZZ,color=col[3])
line25, = axs[1].plot(mean_dS_vmxq,ZZ,color=col[4])
line26, = axs[1].plot(mean_dS_qisf*0,ZZ,color=col[7])
axs[1].set_xlabel('(f) Absolute Salinity anomaly (g/kg)',fontsize=16)
axs[1].set_ylabel('Depth (m)',fontsize=14)
axs[1].xaxis.set_ticks_position('top')
axs[1].xaxis.set_label_position('top')
axs[1].invert_yaxis()
axs[1].set_xlim(-5.2,5.2)
axs[1].set_ylim(1250,0)
axs[1].tick_params(axis='x', labelsize=12)
axs[1].tick_params(axis='y', labelsize=12)
axs[1].grid(True) 
axs[1].text(0.7,700,'Without ice-shelf cavities',color='k',fontsize=14,fontweight='bold',HorizontalAlignment='left')
axs[1].legend((line21,line22,line23,line24,line25,line26),\
              ('$\Delta$ S','Horizontal Advection','Vertical Advection','Horizontal mixing','Vert. Mix. + Surf. Fluxes','Ice Shelf melt'),fontsize=14)

fig.savefig('budget_terms_ContShelf_noisf.jpg')
fig.savefig('budget_terms_ContShelf_noisf.pdf')
