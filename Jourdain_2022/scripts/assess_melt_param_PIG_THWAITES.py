import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os
import functions_nico as nico

#====================================================================

param_type = 'PIGL'  # 'MeanAnt' or 'PIGL'

#====================================================================

def rm_nan(x):
   y = x[~np.isnan(x)]
   return y

#====================================================================

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
ncisfA = ncisfA.rename({'lon':'x'})
ncisfA = ncisfA.rename({'lat':'y'})
ncisfB = xr.open_dataset(file_isfB,decode_cf=False,concat_characters=True)
ncisfB = ncisfB.rename({'lon':'x'})
ncisfB = ncisfB.rename({'lat':'y'})
#ncisfC = xr.open_dataset(file_isfC,decode_cf=False)

isflist = np.array(['Getz', 'Dotson', 'Crosson', 'Thwaites', 'Pine Island', 'Cosgrove' ], dtype=str)
isfnum  = np.array([ 22-1 ,  52-1   ,  38-1    ,  23-1     ,  66-1        ,  53-1      ], dtype=int)

Ncav = np.size(isfnum)
mz = ncmshA.tmask.shape[1]
ZZ=ncmshA.gdept_1d.isel(t=0)

#====================================================================
# load T,S profiles in ice shelf boxes :

npzfile = np.load('TS_profiles_boxes_param.npz')

mean_TpA = npzfile['mean_TpA']
mean_TpB = npzfile['mean_TpB']
mean_TpC = npzfile['mean_TpC']
mean_TfA = npzfile['mean_TfA']
mean_TfB = npzfile['mean_TfB']
mean_TfC = npzfile['mean_TfC']
mean_SpA = npzfile['mean_SpA']
mean_SpB = npzfile['mean_SpB']
mean_SpC = npzfile['mean_SpC']
mean_SfA = npzfile['mean_SfA']
mean_SfB = npzfile['mean_SfB']
mean_SfC = npzfile['mean_SfC']
isflist2 = npzfile['isflist']
isfnum2 = npzfile['isfnum']
ZZ = npzfile['ZZ']

if ( np.size(isflist2) != np.size(isflist) ):
  print('Check consistency between ice shelf lists in this script and in the npz file     >>>>>>> stop')
  exit() 

#====================================================================
# load CMIP5 T,S anomalies in ice shelf boxes :

npzfileCMIP = np.load('TS_profiles_boxes_param_CMIP5anom.npz')

mean_TfA_CMIP = mean_TpA + npzfileCMIP['mean_Tanom']
mean_TfB_CMIP = mean_TpB + npzfileCMIP['mean_Tanom']
mean_TfC_CMIP = mean_TpC + npzfileCMIP['mean_Tanom']
mean_SfA_CMIP = mean_SpA + npzfileCMIP['mean_Sanom']
mean_SfB_CMIP = mean_SpB + npzfileCMIP['mean_Sanom']
mean_SfC_CMIP = mean_SpC + npzfileCMIP['mean_Sanom']

#====================================================================
# Read NEMO's melt rates :

file_SBC_p_A = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MAR/climato_monthly_AMUXL12-GNJ002_BM02MAR_SBC_1989_2009.nc'
file_SBC_f_A = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MARrcp85/climato_monthly_AMUXL12-GNJ002_BM02MARrcp85_SBC_1989_2009.nc'
file_SBC_p_B = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM03MAR/climato_monthly_AMUXL12-GNJ002_BM03MAR_SBC_1989_2009.nc'
file_SBC_f_B = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM03MARrcBDY/climato_monthly_AMUXL12-GNJ002_BM03MARrcBDY_SBC_1989_2009.nc'
file_SBC_p_C = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM04MAR/climato_monthly_AMUXL12-GNJ002_BM04MAR_SBC_1989_2009.nc'
file_SBC_f_C = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM04MARrcp85/climato_monthly_AMUXL12-GNJ002_BM04MARrcp85_SBC_1989_2009.nc'

ncMLTpA = xr.open_dataset(file_SBC_p_A,decode_cf=False)
ncMLTpB = xr.open_dataset(file_SBC_p_B,decode_cf=False)
ncMLTpC = xr.open_dataset(file_SBC_p_C,decode_cf=False)
ncMLTfA = xr.open_dataset(file_SBC_f_A,decode_cf=False)
ncMLTfB = xr.open_dataset(file_SBC_f_B,decode_cf=False)
ncMLTfC = xr.open_dataset(file_SBC_f_C,decode_cf=False)

#====================================================================

def TF_cav_entrance(T,S,Z,maxdepth):
    """
    Remove T,S values deeper than maxdepth and replace them with last vallid value.
    Then calculate Thermal forcing profile.
    Z and maxdepth are positive downward.
    """
    T = np.where( Z <= maxdepth, T , np.nan )
    S = np.where( Z <= maxdepth, S , np.nan )
    ind = np.where(~np.isnan(T))[0]
    last = ind[-1]
    T[last + 1:] = T[last]
    S[last + 1:] = S[last]
    #--
    lbd1 =   -0.0564     # Jourdain et al. (2017)
    lbd2 =    0.0773     # for CT and AS (TEOS10)
    lbd3 =   -7.8633e-8  #
    rhow = 1026.0
    g    =    9.80665
    TF = T - ( lbd1 * S + lbd2 + lbd3 * rhow * g * Z )
    return TF

#====================================================================

# For the KDE distrib of the melt profile:
nKDE=200 # nb of points vertically for the KDE distribution

zKDE=np.zeros((Ncav,nKDE))
distrib_pA_pct05=np.zeros((Ncav,nKDE)) ; distrib_pA_pct50=np.zeros((Ncav,nKDE)) ; distrib_pA_pct95=np.zeros((Ncav,nKDE))
distrib_fA_pct05=np.zeros((Ncav,nKDE)) ; distrib_fA_pct50=np.zeros((Ncav,nKDE)) ; distrib_fA_pct95=np.zeros((Ncav,nKDE))
distrib_pB_pct05=np.zeros((Ncav,nKDE)) ; distrib_pB_pct50=np.zeros((Ncav,nKDE)) ; distrib_pB_pct95=np.zeros((Ncav,nKDE))
distrib_fB_pct05=np.zeros((Ncav,nKDE)) ; distrib_fB_pct50=np.zeros((Ncav,nKDE)) ; distrib_fB_pct95=np.zeros((Ncav,nKDE))
distrib_pC_pct05=np.zeros((Ncav,nKDE)) ; distrib_pC_pct50=np.zeros((Ncav,nKDE)) ; distrib_pC_pct95=np.zeros((Ncav,nKDE))
distrib_fC_pct05=np.zeros((Ncav,nKDE)) ; distrib_fC_pct50=np.zeros((Ncav,nKDE)) ; distrib_fC_pct95=np.zeros((Ncav,nKDE))
NEMOdistrib_pA=np.zeros((Ncav,nKDE))
NEMOdistrib_fA=np.zeros((Ncav,nKDE))
NEMOdistrib_pB=np.zeros((Ncav,nKDE))
NEMOdistrib_fB=np.zeros((Ncav,nKDE))
NEMOdistrib_pC=np.zeros((Ncav,nKDE))
NEMOdistrib_fC=np.zeros((Ncav,nKDE))
distrib_fA_pct50_CMIP=np.zeros((Ncav,nKDE))
distrib_fB_pct50_CMIP=np.zeros((Ncav,nKDE))
distrib_fC_pct50_CMIP=np.zeros((Ncav,nKDE))

rmse_p_pct50=np.zeros((Ncav))
rmse_f_pct50=np.zeros((Ncav))
rmse_f_pct50_CMIP=np.zeros((Ncav))

maxisfdep = np.zeros((Ncav))

# mean square error
def mse(a,b,mask):
  tjp = (a-b)**2 * mask
  mse = tjp.sum() / mask.sum()
  return mse.values

for kcav in np.arange(Ncav):

   print(' ')
   print(isflist[kcav], isfnum[kcav])

   iminA = ncisfA.imin.isel(Nisf=isfnum[kcav]).values.astype('int')-1 ; iminB = ncisfB.imin.isel(Nisf=isfnum[kcav]).values.astype('int')-1 ; iminC = iminB
   imaxA = ncisfA.imax.isel(Nisf=isfnum[kcav]).values.astype('int')-1 ; imaxB = ncisfB.imax.isel(Nisf=isfnum[kcav]).values.astype('int')-1 ; imaxC = imaxB
   jminA = ncisfA.jmin.isel(Nisf=isfnum[kcav]).values.astype('int')-1 ; jminB = ncisfB.jmin.isel(Nisf=isfnum[kcav]).values.astype('int')-1 ; jminC = jminB
   jmaxA = ncisfA.jmax.isel(Nisf=isfnum[kcav]).values.astype('int')-1 ; jmaxB = ncisfB.jmax.isel(Nisf=isfnum[kcav]).values.astype('int')-1 ; jmaxC = jmaxB
   imin = np.min([iminA,iminB,iminC])
   imax = np.max([imaxA,imaxB,imaxC])
   jmin = np.min([jminA,jminB,jminC])
   jmax = np.max([jmaxA,jmaxB,jmaxC])

   mskA = ncmshA.tmask.max(axis=1).isel(t=0,x=slice(imin,imax),y=slice(jmin,jmax)).where( ( ncisfA.isfmask.isel(x=slice(imin,imax),y=slice(jmin,jmax)) == isfnum[kcav]+1 ) , 0.e0 )
   mskB = ncmshB.tmask.max(axis=1).isel(t=0,x=slice(imin,imax),y=slice(jmin,jmax)).where( ( ncisfB.isfmask.isel(x=slice(imin,imax),y=slice(jmin,jmax)) == isfnum[kcav]+1 ) , 0.e0 )
   mskC = mskB

   ZdA = ncmshA.isfdraft.isel(t=0,x=slice(imin,imax),y=slice(jmin,jmax)) * mskA
   ZdB = ncmshB.isfdraft.isel(t=0,x=slice(imin,imax),y=slice(jmin,jmax)) * mskB
   ZdC = ZdB
   maxisfdep[kcav] = np.max( [ ZdA.max().values, ZdB.max().values ] )
   print('max ice draft depth : ', maxisfdep[kcav], ' m')

   areacella = ncmshA.e1t.isel(t=0,x=slice(imin,imax),y=slice(jmin,jmax)) * ncmshA.e2t.isel(t=0,x=slice(imin,imax),y=slice(jmin,jmax))

   print('(imin, imax, jmin, jmax) = ', imin, imax, jmin, jmax)
   #plt.pcolormesh(ncmshB.misf.isel(t=0,x=slice(imin,imax),y=slice(jmin,jmax)))
   #plt.pcolormesh(ncisfB.isfmask.isel(x=slice(imin,imax),y=slice(jmin,jmax)))
   #plt.pcolormesh(mskA)
   #plt.colorbar()
   #plt.show()

   #--------------------------------------------
   # melt simulated by NEMO:
   NEMOfac = -86400. * 365.25 * 1.e-3 # kg /m2 / s  ->  m.w.e /yr
   NEMOmelt_pA = ncMLTpA.fwfisf.mean(axis=0).isel(x=slice(imin,imax),y=slice(jmin,jmax)) * mskA * NEMOfac
   NEMOmelt_fA = ncMLTfA.fwfisf.mean(axis=0).isel(x=slice(imin,imax),y=slice(jmin,jmax)) * mskA * NEMOfac
   NEMOmelt_pB = ncMLTpB.fwfisf.mean(axis=0).isel(x=slice(imin,imax),y=slice(jmin,jmax)) * mskB * NEMOfac
   NEMOmelt_fB = ncMLTfB.fwfisf.mean(axis=0).isel(x=slice(imin,imax),y=slice(jmin,jmax)) * mskB * NEMOfac
   NEMOmelt_pC = ncMLTpC.fwfisf.mean(axis=0).isel(x=slice(imin,imax),y=slice(jmin,jmax)) * mskC * NEMOfac
   NEMOmelt_fC = ncMLTfC.fwfisf.mean(axis=0).isel(x=slice(imin,imax),y=slice(jmin,jmax)) * mskC * NEMOfac

   #--------------------------------------------------------
   # Vertical profiles of thermal forcing from NEMO's T,S:
   # we don't keep input T,S for depth deeper then the maximum seafloor depth of the entrance :
   TFpA = TF_cav_entrance( mean_TpA[kcav,:], mean_SpA[kcav,:], ZZ, ncisfA.front_bot_dep_max.isel(Nisf=isfnum[kcav]).values )
   TFpB = TF_cav_entrance( mean_TpB[kcav,:], mean_SpB[kcav,:], ZZ, ncisfB.front_bot_dep_max.isel(Nisf=isfnum[kcav]).values )
   TFpC = TF_cav_entrance( mean_TpC[kcav,:], mean_SpC[kcav,:], ZZ, ncisfB.front_bot_dep_max.isel(Nisf=isfnum[kcav]).values )
   #--
   TFfA = TF_cav_entrance( mean_TfA[kcav,:], mean_SfA[kcav,:], ZZ, ncisfA.front_bot_dep_max.isel(Nisf=isfnum[kcav]).values )
   TFfB = TF_cav_entrance( mean_TfB[kcav,:], mean_SfB[kcav,:], ZZ, ncisfB.front_bot_dep_max.isel(Nisf=isfnum[kcav]).values )
   TFfC = TF_cav_entrance( mean_TfC[kcav,:], mean_SfC[kcav,:], ZZ, ncisfB.front_bot_dep_max.isel(Nisf=isfnum[kcav]).values )
   #--
   TFdraft_pA = np.interp( ZdA, ZZ, TFpA )
   TFdraft_fA = np.interp( ZdA, ZZ, TFfA )
   TFdraft_pB = np.interp( ZdB, ZZ, TFpB )
   TFdraft_fB = np.interp( ZdB, ZZ, TFfB )
   TFdraft_pC = np.interp( ZdC, ZZ, TFpC )
   TFdraft_fC = np.interp( ZdC, ZZ, TFfC )
   #--
   mean_TFdraft_pA = np.sum( np.sum( TFdraft_pA * areacella.values * mskA.values , axis=1 ), axis=0 ) / np.sum( np.sum( areacella.values * mskA.values , axis=1 ), axis=0 )
   mean_TFdraft_fA = np.sum( np.sum( TFdraft_fA * areacella.values * mskA.values , axis=1 ), axis=0 ) / np.sum( np.sum( areacella.values * mskA.values , axis=1 ), axis=0 )
   mean_TFdraft_pB = np.sum( np.sum( TFdraft_pB * areacella.values * mskB.values , axis=1 ), axis=0 ) / np.sum( np.sum( areacella.values * mskB.values , axis=1 ), axis=0 )
   mean_TFdraft_fB = np.sum( np.sum( TFdraft_fB * areacella.values * mskB.values , axis=1 ), axis=0 ) / np.sum( np.sum( areacella.values * mskB.values , axis=1 ), axis=0 )
   mean_TFdraft_pC = np.sum( np.sum( TFdraft_pC * areacella.values * mskC.values , axis=1 ), axis=0 ) / np.sum( np.sum( areacella.values * mskC.values , axis=1 ), axis=0 )
   mean_TFdraft_fC = np.sum( np.sum( TFdraft_fC * areacella.values * mskC.values , axis=1 ), axis=0 ) / np.sum( np.sum( areacella.values * mskC.values , axis=1 ), axis=0 )
   print('TF NEMO: ', mean_TFdraft_pA, mean_TFdraft_fA, mean_TFdraft_pB, mean_TFdraft_fB, mean_TFdraft_pC, mean_TFdraft_fC )

   #----------------------------------------------------------------
   # Vertical profiles of thermal forcing from NEMO + CMIP5 anomaly :
   # we don't keep input T,S for depth deeper then the maximum seafloor depth of the entrance :
   TFfA_CMIP = TF_cav_entrance( mean_TfA_CMIP[kcav,:], mean_SfA_CMIP[kcav,:], ZZ, ncisfA.front_bot_dep_max.isel(Nisf=isfnum[kcav]).values )
   TFfB_CMIP = TF_cav_entrance( mean_TfB_CMIP[kcav,:], mean_SfB_CMIP[kcav,:], ZZ, ncisfB.front_bot_dep_max.isel(Nisf=isfnum[kcav]).values )
   TFfC_CMIP = TF_cav_entrance( mean_TfC_CMIP[kcav,:], mean_SfC_CMIP[kcav,:], ZZ, ncisfB.front_bot_dep_max.isel(Nisf=isfnum[kcav]).values )
   #--
   TFdraft_fA_CMIP = np.interp( ZdA, ZZ, TFfA_CMIP )
   TFdraft_fB_CMIP = np.interp( ZdB, ZZ, TFfB_CMIP )
   TFdraft_fC_CMIP = np.interp( ZdC, ZZ, TFfC_CMIP )
   #--
   mean_TFdraft_fA_CMIP = np.sum( np.sum( TFdraft_fA_CMIP * areacella.values * mskA.values , axis=1 ), axis=0 ) / np.sum( np.sum( areacella.values * mskA.values , axis=1 ), axis=0 )
   mean_TFdraft_fB_CMIP = np.sum( np.sum( TFdraft_fB_CMIP * areacella.values * mskB.values , axis=1 ), axis=0 ) / np.sum( np.sum( areacella.values * mskB.values , axis=1 ), axis=0 )
   mean_TFdraft_fC_CMIP = np.sum( np.sum( TFdraft_fC_CMIP * areacella.values * mskC.values , axis=1 ), axis=0 ) / np.sum( np.sum( areacella.values * mskC.values , axis=1 ), axis=0 )
   print('TF CMIP anom: ', mean_TFdraft_fA_CMIP, mean_TFdraft_fB_CMIP, mean_TFdraft_fC_CMIP )

   #---------
   # Quadratic semi-local (ISMIP6):
   if ( param_type == 'MeanAnt' ):
     gamma0_pct05 =  9620.
     gamma0_pct50 = 14500.
     gamma0_pct95 = 21000.
     deltaT = 0.0 # empirical correction
   elif ( param_type == 'PIGL' ):
     gamma0_pct05 =  88000.
     gamma0_pct50 = 159000.
     gamma0_pct95 = 471000.
     deltaT = -1.90 # empirical correction
   K2 = (1028. * 3974. / ( 918. * 3.34e5 ) )**2 # ((rhow.cpw)/(rhoi.Lf))^2
   #--
   melt_pA_pct05 = gamma0_pct05 * K2 * ( TFdraft_pA + deltaT ) * np.abs( mean_TFdraft_pA + deltaT ) 
   melt_fA_pct05 = gamma0_pct05 * K2 * ( TFdraft_fA + deltaT ) * np.abs( mean_TFdraft_fA + deltaT ) 
   melt_pB_pct05 = gamma0_pct05 * K2 * ( TFdraft_pB + deltaT ) * np.abs( mean_TFdraft_pB + deltaT ) 
   melt_fB_pct05 = gamma0_pct05 * K2 * ( TFdraft_fB + deltaT ) * np.abs( mean_TFdraft_fB + deltaT ) 
   melt_pC_pct05 = gamma0_pct05 * K2 * ( TFdraft_pC + deltaT ) * np.abs( mean_TFdraft_pC + deltaT ) 
   melt_fC_pct05 = gamma0_pct05 * K2 * ( TFdraft_fC + deltaT ) * np.abs( mean_TFdraft_fC + deltaT )
   #--
   melt_pA_pct50 = gamma0_pct50 * K2 * ( TFdraft_pA + deltaT ) * np.abs( mean_TFdraft_pA + deltaT )
   melt_fA_pct50 = gamma0_pct50 * K2 * ( TFdraft_fA + deltaT ) * np.abs( mean_TFdraft_fA + deltaT )
   melt_pB_pct50 = gamma0_pct50 * K2 * ( TFdraft_pB + deltaT ) * np.abs( mean_TFdraft_pB + deltaT )
   melt_fB_pct50 = gamma0_pct50 * K2 * ( TFdraft_fB + deltaT ) * np.abs( mean_TFdraft_fB + deltaT )
   melt_pC_pct50 = gamma0_pct50 * K2 * ( TFdraft_pC + deltaT ) * np.abs( mean_TFdraft_pC + deltaT )
   melt_fC_pct50 = gamma0_pct50 * K2 * ( TFdraft_fC + deltaT ) * np.abs( mean_TFdraft_fC + deltaT )
   #--
   melt_pA_pct95 = gamma0_pct95 * K2 * ( TFdraft_pA + deltaT ) * np.abs( mean_TFdraft_pA + deltaT )
   melt_fA_pct95 = gamma0_pct95 * K2 * ( TFdraft_fA + deltaT ) * np.abs( mean_TFdraft_fA + deltaT )
   melt_pB_pct95 = gamma0_pct95 * K2 * ( TFdraft_pB + deltaT ) * np.abs( mean_TFdraft_pB + deltaT )
   melt_fB_pct95 = gamma0_pct95 * K2 * ( TFdraft_fB + deltaT ) * np.abs( mean_TFdraft_fB + deltaT )
   melt_pC_pct95 = gamma0_pct95 * K2 * ( TFdraft_pC + deltaT ) * np.abs( mean_TFdraft_pC + deltaT )
   melt_fC_pct95 = gamma0_pct95 * K2 * ( TFdraft_fC + deltaT ) * np.abs( mean_TFdraft_fC + deltaT )
   #--
   melt_fA_pct50_CMIP = gamma0_pct50 * K2 * ( TFdraft_fA_CMIP + deltaT ) * np.abs( mean_TFdraft_fA_CMIP + deltaT )
   melt_fB_pct50_CMIP = gamma0_pct50 * K2 * ( TFdraft_fB_CMIP + deltaT ) * np.abs( mean_TFdraft_fB_CMIP + deltaT )
   melt_fC_pct50_CMIP = gamma0_pct50 * K2 * ( TFdraft_fC_CMIP + deltaT ) * np.abs( mean_TFdraft_fC_CMIP + deltaT )

   #----------
   # RMSE

   rmse_p_pct50[kcav] = np.sqrt( mse(melt_pA_pct50,NEMOmelt_pA,mskA) + mse(melt_pB_pct50,NEMOmelt_pB,mskB) + mse(melt_pC_pct50,NEMOmelt_pC,mskC) ) / np.sqrt(3.)
   rmse_f_pct50[kcav] = np.sqrt( mse(melt_fA_pct50,NEMOmelt_fA,mskA) + mse(melt_fB_pct50,NEMOmelt_fB,mskB) + mse(melt_fC_pct50,NEMOmelt_fC,mskC) ) / np.sqrt(3.)

   rmse_f_pct50_CMIP[kcav] = np.sqrt( mse(melt_fA_pct50_CMIP,NEMOmelt_fA,mskA) + mse(melt_fB_pct50_CMIP,NEMOmelt_fB,mskB) + mse(melt_fC_pct50_CMIP,NEMOmelt_fC,mskC) ) / np.sqrt(3.)

   print('RMSE : ',rmse_p_pct50[kcav], rmse_f_pct50[kcav], rmse_f_pct50_CMIP[kcav] )
 
   #----------
   # KDE fit (melt distribution) of vertical melt profiles :

   #--
   melt_pA_vec_pct05 = np.where( (mskA.values>0), melt_pA_pct05, np.nan ) ; melt_pA_vec_pct05 = rm_nan(melt_pA_vec_pct05)
   melt_fA_vec_pct05 = np.where( (mskA.values>0), melt_fA_pct05, np.nan ) ; melt_fA_vec_pct05 = rm_nan(melt_fA_vec_pct05)
   melt_pB_vec_pct05 = np.where( (mskB.values>0), melt_pB_pct05, np.nan ) ; melt_pB_vec_pct05 = rm_nan(melt_pB_vec_pct05)
   melt_fB_vec_pct05 = np.where( (mskB.values>0), melt_fB_pct05, np.nan ) ; melt_fB_vec_pct05 = rm_nan(melt_fB_vec_pct05)
   melt_pC_vec_pct05 = np.where( (mskC.values>0), melt_pC_pct05, np.nan ) ; melt_pC_vec_pct05 = rm_nan(melt_pC_vec_pct05)
   melt_fC_vec_pct05 = np.where( (mskC.values>0), melt_fC_pct05, np.nan ) ; melt_fC_vec_pct05 = rm_nan(melt_fC_vec_pct05)
   #--
   melt_pA_vec_pct50 = np.where( (mskA.values>0), melt_pA_pct50, np.nan ) ; melt_pA_vec_pct50 = rm_nan(melt_pA_vec_pct50)
   melt_fA_vec_pct50 = np.where( (mskA.values>0), melt_fA_pct50, np.nan ) ; melt_fA_vec_pct50 = rm_nan(melt_fA_vec_pct50)
   melt_pB_vec_pct50 = np.where( (mskB.values>0), melt_pB_pct50, np.nan ) ; melt_pB_vec_pct50 = rm_nan(melt_pB_vec_pct50)
   melt_fB_vec_pct50 = np.where( (mskB.values>0), melt_fB_pct50, np.nan ) ; melt_fB_vec_pct50 = rm_nan(melt_fB_vec_pct50)
   melt_pC_vec_pct50 = np.where( (mskC.values>0), melt_pC_pct50, np.nan ) ; melt_pC_vec_pct50 = rm_nan(melt_pC_vec_pct50)
   melt_fC_vec_pct50 = np.where( (mskC.values>0), melt_fC_pct50, np.nan ) ; melt_fC_vec_pct50 = rm_nan(melt_fC_vec_pct50)
   #--
   melt_pA_vec_pct95 = np.where( (mskA.values>0), melt_pA_pct95, np.nan ) ; melt_pA_vec_pct95 = rm_nan(melt_pA_vec_pct95)
   melt_fA_vec_pct95 = np.where( (mskA.values>0), melt_fA_pct95, np.nan ) ; melt_fA_vec_pct95 = rm_nan(melt_fA_vec_pct95)
   melt_pB_vec_pct95 = np.where( (mskB.values>0), melt_pB_pct95, np.nan ) ; melt_pB_vec_pct95 = rm_nan(melt_pB_vec_pct95)
   melt_fB_vec_pct95 = np.where( (mskB.values>0), melt_fB_pct95, np.nan ) ; melt_fB_vec_pct95 = rm_nan(melt_fB_vec_pct95)
   melt_pC_vec_pct95 = np.where( (mskC.values>0), melt_pC_pct95, np.nan ) ; melt_pC_vec_pct95 = rm_nan(melt_pC_vec_pct95)
   melt_fC_vec_pct95 = np.where( (mskC.values>0), melt_fC_pct95, np.nan ) ; melt_fC_vec_pct95 = rm_nan(melt_fC_vec_pct95)
   #--
   melt_fA_vec_pct50_CMIP = np.where( (mskA.values>0), melt_fA_pct50_CMIP, np.nan ) ; melt_fA_vec_pct50_CMIP = rm_nan(melt_fA_vec_pct50_CMIP)
   melt_fB_vec_pct50_CMIP = np.where( (mskB.values>0), melt_fB_pct50_CMIP, np.nan ) ; melt_fB_vec_pct50_CMIP = rm_nan(melt_fB_vec_pct50_CMIP)
   melt_fC_vec_pct50_CMIP = np.where( (mskC.values>0), melt_fC_pct50_CMIP, np.nan ) ; melt_fC_vec_pct50_CMIP = rm_nan(melt_fC_vec_pct50_CMIP)
   #--
   NEMOmelt_pA_vec = np.where( (mskA.values>0), NEMOmelt_pA, np.nan ) ; NEMOmelt_pA_vec = rm_nan(NEMOmelt_pA_vec)
   NEMOmelt_fA_vec = np.where( (mskA.values>0), NEMOmelt_fA, np.nan ) ; NEMOmelt_fA_vec = rm_nan(NEMOmelt_fA_vec)
   NEMOmelt_pB_vec = np.where( (mskB.values>0), NEMOmelt_pB, np.nan ) ; NEMOmelt_pB_vec = rm_nan(NEMOmelt_pB_vec)
   NEMOmelt_fB_vec = np.where( (mskB.values>0), NEMOmelt_fB, np.nan ) ; NEMOmelt_fB_vec = rm_nan(NEMOmelt_fB_vec)
   NEMOmelt_pC_vec = np.where( (mskC.values>0), NEMOmelt_pC, np.nan ) ; NEMOmelt_pC_vec = rm_nan(NEMOmelt_pC_vec)
   NEMOmelt_fC_vec = np.where( (mskC.values>0), NEMOmelt_fC, np.nan ) ; NEMOmelt_fC_vec = rm_nan(NEMOmelt_fC_vec)
   #--
   ZdA_vec = np.where( (mskA.values>0), ZdA, np.nan ) ; ZdA_vec = rm_nan(ZdA_vec)
   ZdB_vec = np.where( (mskB.values>0), ZdB, np.nan ) ; ZdB_vec = rm_nan(ZdB_vec)
   ZdC_vec = np.where( (mskC.values>0), ZdC, np.nan ) ; ZdC_vec = rm_nan(ZdC_vec)
   print( np.amax(ZdA_vec), np.amax(ZdB_vec), np.amax(ZdC_vec) )
   #--
   zKDE[kcav,:] = np.linspace(30.,maxisfdep[kcav],nKDE)
   zsigma=maxisfdep[kcav]/20
   for kk in np.arange(0,nKDE):
     kernelA = np.exp(-(ZdA_vec - zKDE[kcav,kk]) ** 2 / (2 * zsigma ** 2)) ;  kernelA = kernelA / sum(kernelA)
     kernelB = np.exp(-(ZdB_vec - zKDE[kcav,kk]) ** 2 / (2 * zsigma ** 2)) ;  kernelB = kernelB / sum(kernelB)
     kernelC = np.exp(-(ZdC_vec - zKDE[kcav,kk]) ** 2 / (2 * zsigma ** 2)) ;  kernelC = kernelC / sum(kernelC)
     #--
     distrib_pA_pct05[kcav,kk] = sum( melt_pA_vec_pct05 * kernelA )
     distrib_fA_pct05[kcav,kk] = sum( melt_fA_vec_pct05 * kernelA )
     distrib_pB_pct05[kcav,kk] = sum( melt_pB_vec_pct05 * kernelB )
     distrib_fB_pct05[kcav,kk] = sum( melt_fB_vec_pct05 * kernelB )
     distrib_pC_pct05[kcav,kk] = sum( melt_pC_vec_pct05 * kernelC )
     distrib_fC_pct05[kcav,kk] = sum( melt_fC_vec_pct05 * kernelC )
     #--
     distrib_pA_pct50[kcav,kk] = sum( melt_pA_vec_pct50 * kernelA )
     distrib_fA_pct50[kcav,kk] = sum( melt_fA_vec_pct50 * kernelA )
     distrib_pB_pct50[kcav,kk] = sum( melt_pB_vec_pct50 * kernelB )
     distrib_fB_pct50[kcav,kk] = sum( melt_fB_vec_pct50 * kernelB )
     distrib_pC_pct50[kcav,kk] = sum( melt_pC_vec_pct50 * kernelC )
     distrib_fC_pct50[kcav,kk] = sum( melt_fC_vec_pct50 * kernelC )
     #--
     distrib_pA_pct95[kcav,kk] = sum( melt_pA_vec_pct95 * kernelA )
     distrib_fA_pct95[kcav,kk] = sum( melt_fA_vec_pct95 * kernelA )
     distrib_pB_pct95[kcav,kk] = sum( melt_pB_vec_pct95 * kernelB )
     distrib_fB_pct95[kcav,kk] = sum( melt_fB_vec_pct95 * kernelB )
     distrib_pC_pct95[kcav,kk] = sum( melt_pC_vec_pct95 * kernelC )
     distrib_fC_pct95[kcav,kk] = sum( melt_fC_vec_pct95 * kernelC )
     #--
     NEMOdistrib_pA[kcav,kk] = sum( NEMOmelt_pA_vec * kernelA )
     NEMOdistrib_fA[kcav,kk] = sum( NEMOmelt_fA_vec * kernelA )
     NEMOdistrib_pB[kcav,kk] = sum( NEMOmelt_pB_vec * kernelB )
     NEMOdistrib_fB[kcav,kk] = sum( NEMOmelt_fB_vec * kernelB )
     NEMOdistrib_pC[kcav,kk] = sum( NEMOmelt_pC_vec * kernelC )
     NEMOdistrib_fC[kcav,kk] = sum( NEMOmelt_fC_vec * kernelC )
     #--
     distrib_fA_pct50_CMIP[kcav,kk] = sum( melt_fA_vec_pct50_CMIP * kernelA )
     distrib_fB_pct50_CMIP[kcav,kk] = sum( melt_fB_vec_pct50_CMIP * kernelB )
     distrib_fC_pct50_CMIP[kcav,kk] = sum( melt_fC_vec_pct50_CMIP * kernelC )

#================================================================
# Plots :

fig, axs = plt.subplots(nrows=1,ncols=2,figsize=(21.0,9.9))
axs = axs.ravel()

#----------
# Pine Island:
kcav=4
shadpA, = axs[0].fill( np.concatenate( [distrib_pA_pct05[kcav,:],np.flipud(distrib_pA_pct95[kcav,:])] ),\
                       np.concatenate( [zKDE[kcav,:],np.flipud(zKDE[kcav,:])] ),\
                       'cornflowerblue',alpha=0.1 )
shadfA, = axs[0].fill( np.concatenate( [distrib_fA_pct05[kcav,:],np.flipud(distrib_fA_pct95[kcav,:])] ),\
                       np.concatenate( [zKDE[kcav,:],np.flipud(zKDE[kcav,:])] ),\
                       'orange',alpha=0.1 )
shadpB, = axs[0].fill( np.concatenate( [distrib_pB_pct05[kcav,:],np.flipud(distrib_pB_pct95[kcav,:])] ),\
                       np.concatenate( [zKDE[kcav,:],np.flipud(zKDE[kcav,:])] ),\
                       'cornflowerblue',alpha=0.1 )
shadfB, = axs[0].fill( np.concatenate( [distrib_fB_pct05[kcav,:],np.flipud(distrib_fB_pct95[kcav,:])] ),\
                       np.concatenate( [zKDE[kcav,:],np.flipud(zKDE[kcav,:])] ),\
                       'orange',alpha=0.1 )
shadpC, = axs[0].fill( np.concatenate( [distrib_pC_pct05[kcav,:],np.flipud(distrib_pC_pct95[kcav,:])] ),\
                       np.concatenate( [zKDE[kcav,:],np.flipud(zKDE[kcav,:])] ),\
                       'cornflowerblue',alpha=0.1 )
shadfC, = axs[0].fill( np.concatenate( [distrib_fC_pct05[kcav,:],np.flipud(distrib_fC_pct95[kcav,:])] ),\
                       np.concatenate( [zKDE[kcav,:],np.flipud(zKDE[kcav,:])] ),\
                       'orange',alpha=0.1 )

linepA, = axs[0].plot(distrib_pA_pct50[kcav,:],zKDE[kcav,:],color='cornflowerblue',linestyle=(10,(4, 6)),linewidth=1.2)
linefA, = axs[0].plot(distrib_fA_pct50[kcav,:],zKDE[kcav,:],color='orange',linestyle=(10,(4, 6)),linewidth=1.2)
linepB, = axs[0].plot(distrib_pB_pct50[kcav,:],zKDE[kcav,:],color='cornflowerblue',linestyle=(10,(4, 6)),linewidth=1.2)
linefB, = axs[0].plot(distrib_fB_pct50[kcav,:],zKDE[kcav,:],color='orange',linestyle=(10,(4, 6)),linewidth=1.2)
linepC, = axs[0].plot(distrib_pC_pct50[kcav,:],zKDE[kcav,:],color='cornflowerblue',linestyle=(10,(4, 6)),linewidth=1.2)
linefC, = axs[0].plot(distrib_fC_pct50[kcav,:],zKDE[kcav,:],color='orange',linestyle=(10,(4, 6)),linewidth=1.2)

lineNEMOpA, = axs[0].plot(NEMOdistrib_pA[kcav,:],zKDE[kcav,:],color='cornflowerblue',linewidth=1.5)
lineNEMOpB, = axs[0].plot(NEMOdistrib_pB[kcav,:],zKDE[kcav,:],color='cornflowerblue',linewidth=1.5)
lineNEMOpC, = axs[0].plot(NEMOdistrib_pC[kcav,:],zKDE[kcav,:],color='cornflowerblue',linewidth=1.5)
lineNEMOfA, = axs[0].plot(NEMOdistrib_fA[kcav,:],zKDE[kcav,:],color='orange',linewidth=1.5)
lineNEMOfB, = axs[0].plot(NEMOdistrib_fB[kcav,:],zKDE[kcav,:],color='orange',linewidth=1.5)
lineNEMOfC, = axs[0].plot(NEMOdistrib_fC[kcav,:],zKDE[kcav,:],color='orange',linewidth=1.5)

lineCMIPfA, = axs[0].plot(distrib_fA_pct50_CMIP[kcav,:],zKDE[kcav,:],color='brown',linestyle=(10,(4, 6)),linewidth=1.2)
lineCMIPfB, = axs[0].plot(distrib_fB_pct50_CMIP[kcav,:],zKDE[kcav,:],color='brown',linestyle=(10,(4, 6)),linewidth=1.2)
lineCMIPfC, = axs[0].plot(distrib_fC_pct50_CMIP[kcav,:],zKDE[kcav,:],color='brown',linestyle=(10,(4, 6)),linewidth=1.2)

if ( param_type == 'MeanAnt' ):
  lab = '(a) '+isflist[kcav]+' ('+param_type+' calibration)'
else:
  lab = '(c) '+isflist[kcav]+' ('+param_type+' calibration)'
axs[0].set_title(lab,fontsize=24,weight='bold')
axs[0].set_xlabel('melt rate (m/yr)',fontsize=20)
axs[0].set_ylabel('Depth (m)',fontsize=20)
axs[0].xaxis.set_ticks_position('top')
axs[0].xaxis.set_label_position('top')
axs[0].invert_yaxis()
if ( param_type == 'MeanAnt' ):
  axs[0].set_xlim(left=0.)
axs[0].set_ylim(maxisfdep[kcav],30)
axs[0].tick_params(axis='x', labelsize=16)
axs[0].tick_params(axis='y', labelsize=16)
if ( param_type == 'MeanAnt' ):
  axs[0].legend((lineNEMOpA,lineNEMOfA,linepA,shadpA,linefA,shadfA,lineCMIPfA),\
                ('NEMO\'s present-day','NEMO\'s future','param. present-day','5th - 95th range','param. future','5th - 95th range','param. fut. (CMIP5 anom.)'),fontsize=17)

[xmin, xmax]=axs[0].get_xlim()
[ymin, ymax]=axs[0].get_ylim()
if ( param_type == 'MeanAnt' ):
  string='RMSE = '+np.array2string(rmse_p_pct50[kcav],precision=1,floatmode='fixed')
  axs[0].text(xmin+0.70*(xmax-xmin),ymin+0.23*(ymax-ymin),string,color='cornflowerblue',fontsize=21)
  string='RMSE = '+np.array2string(rmse_f_pct50[kcav],precision=1,floatmode='fixed')
  axs[0].text(xmin+0.70*(xmax-xmin),ymin+0.15*(ymax-ymin),string,color='orange',fontsize=21)
  string='RMSE = '+np.array2string(rmse_f_pct50_CMIP[kcav],precision=1,floatmode='fixed')
  axs[0].text(xmin+0.70*(xmax-xmin),ymin+0.07*(ymax-ymin),string,color='brown',fontsize=21)
else:
  string='RMSE = '+np.array2string(rmse_p_pct50[kcav],precision=1,floatmode='fixed')
  axs[0].text(xmin+0.70*(xmax-xmin),ymin+0.93*(ymax-ymin),string,color='cornflowerblue',fontsize=21)
  string='RMSE = '+np.array2string(rmse_f_pct50[kcav],precision=1,floatmode='fixed')
  axs[0].text(xmin+0.70*(xmax-xmin),ymin+0.85*(ymax-ymin),string,color='orange',fontsize=21)
  string='RMSE = '+np.array2string(rmse_f_pct50_CMIP[kcav],precision=1,floatmode='fixed')
  axs[0].text(xmin+0.70*(xmax-xmin),ymin+0.77*(ymax-ymin),string,color='brown',fontsize=21)

#----------
# Thwaites :
kcav=3
shadpA, = axs[1].fill( np.concatenate( [distrib_pA_pct05[kcav,:],np.flipud(distrib_pA_pct95[kcav,:])] ),\
                       np.concatenate( [zKDE[kcav,:],np.flipud(zKDE[kcav,:])] ),\
                       'cornflowerblue',alpha=0.1 )
shadfA, = axs[1].fill( np.concatenate( [distrib_fA_pct05[kcav,:],np.flipud(distrib_fA_pct95[kcav,:])] ),\
                       np.concatenate( [zKDE[kcav,:],np.flipud(zKDE[kcav,:])] ),\
                       'orange',alpha=0.1 )
shadpB, = axs[1].fill( np.concatenate( [distrib_pB_pct05[kcav,:],np.flipud(distrib_pB_pct95[kcav,:])] ),\
                       np.concatenate( [zKDE[kcav,:],np.flipud(zKDE[kcav,:])] ),\
                       'cornflowerblue',alpha=0.1 )
shadfB, = axs[1].fill( np.concatenate( [distrib_fB_pct05[kcav,:],np.flipud(distrib_fB_pct95[kcav,:])] ),\
                       np.concatenate( [zKDE[kcav,:],np.flipud(zKDE[kcav,:])] ),\
                       'orange',alpha=0.1 )
shadpC, = axs[1].fill( np.concatenate( [distrib_pC_pct05[kcav,:],np.flipud(distrib_pC_pct95[kcav,:])] ),\
                       np.concatenate( [zKDE[kcav,:],np.flipud(zKDE[kcav,:])] ),\
                       'cornflowerblue',alpha=0.1 )
shadfC, = axs[1].fill( np.concatenate( [distrib_fC_pct05[kcav,:],np.flipud(distrib_fC_pct95[kcav,:])] ),\
                       np.concatenate( [zKDE[kcav,:],np.flipud(zKDE[kcav,:])] ),\
                       'orange',alpha=0.1 )

linepA, = axs[1].plot(distrib_pA_pct50[kcav,:],zKDE[kcav,:],color='cornflowerblue',linestyle=(10,(4, 6)),linewidth=1.2)
linefA, = axs[1].plot(distrib_fA_pct50[kcav,:],zKDE[kcav,:],color='orange',linestyle=(10,(4, 6)),linewidth=1.2)
linepB, = axs[1].plot(distrib_pB_pct50[kcav,:],zKDE[kcav,:],color='cornflowerblue',linestyle=(10,(4, 6)),linewidth=1.2)
linefB, = axs[1].plot(distrib_fB_pct50[kcav,:],zKDE[kcav,:],color='orange',linestyle=(10,(4, 6)),linewidth=1.2)
linepC, = axs[1].plot(distrib_pC_pct50[kcav,:],zKDE[kcav,:],color='cornflowerblue',linestyle=(10,(4, 6)),linewidth=1.2)
linefC, = axs[1].plot(distrib_fC_pct50[kcav,:],zKDE[kcav,:],color='orange',linestyle=(10,(4, 6)),linewidth=1.2)

lineNEMOpA, = axs[1].plot(NEMOdistrib_pA[kcav,:],zKDE[kcav,:],color='cornflowerblue',linewidth=1.5)
lineNEMOpB, = axs[1].plot(NEMOdistrib_pB[kcav,:],zKDE[kcav,:],color='cornflowerblue',linewidth=1.5)
lineNEMOpC, = axs[1].plot(NEMOdistrib_pC[kcav,:],zKDE[kcav,:],color='cornflowerblue',linewidth=1.5)
lineNEMOfA, = axs[1].plot(NEMOdistrib_fA[kcav,:],zKDE[kcav,:],color='orange',linewidth=1.5)
lineNEMOfB, = axs[1].plot(NEMOdistrib_fB[kcav,:],zKDE[kcav,:],color='orange',linewidth=1.5)
lineNEMOfC, = axs[1].plot(NEMOdistrib_fC[kcav,:],zKDE[kcav,:],color='orange',linewidth=1.5)

lineCMIPfA, = axs[1].plot(distrib_fA_pct50_CMIP[kcav,:],zKDE[kcav,:],color='brown',linestyle=(10,(4, 6)),linewidth=1.2)
lineCMIPfB, = axs[1].plot(distrib_fB_pct50_CMIP[kcav,:],zKDE[kcav,:],color='brown',linestyle=(10,(4, 6)),linewidth=1.2)
lineCMIPfC, = axs[1].plot(distrib_fC_pct50_CMIP[kcav,:],zKDE[kcav,:],color='brown',linestyle=(10,(4, 6)),linewidth=1.2)

if ( param_type == 'MeanAnt' ):
  lab = '(b) '+isflist[kcav]+' ('+param_type+' calibration)'
else:
  lab = '(d) '+isflist[kcav]+' ('+param_type+' calibration)'
axs[1].set_title(lab,fontsize=24,weight='bold')
axs[1].set_xlabel('melt rate (m/yr)',fontsize=20)
axs[1].set_ylabel('Depth (m)',fontsize=20)
axs[1].xaxis.set_ticks_position('top')
axs[1].xaxis.set_label_position('top')
axs[1].invert_yaxis()
if ( param_type == 'MeanAnt' ):
  axs[1].set_xlim(left=0.)
axs[1].set_ylim(maxisfdep[kcav],30)
axs[1].tick_params(axis='x', labelsize=16)
axs[1].tick_params(axis='y', labelsize=16)

[xmin, xmax]=axs[1].get_xlim()
[ymin, ymax]=axs[1].get_ylim()
string='RMSE = '+np.array2string(rmse_p_pct50[kcav],precision=1,floatmode='fixed')
axs[1].text(xmin+0.70*(xmax-xmin),ymin+0.93*(ymax-ymin),string,color='cornflowerblue',fontsize=21)
string='RMSE = '+np.array2string(rmse_f_pct50[kcav],precision=1,floatmode='fixed')
axs[1].text(xmin+0.70*(xmax-xmin),ymin+0.85*(ymax-ymin),string,color='orange',fontsize=21)
string='RMSE = '+np.array2string(rmse_f_pct50_CMIP[kcav],precision=1,floatmode='fixed')
axs[1].text(xmin+0.70*(xmax-xmin),ymin+0.77*(ymax-ymin),string,color='brown',fontsize=21)


#----------------
filejpg='melt_param_assessment_PIG_THW_'+param_type+'.jpg'
filepdf='melt_param_assessment_PIG_THW_'+param_type+'.pdf'
fig.savefig(filejpg)
fig.savefig(filepdf)

