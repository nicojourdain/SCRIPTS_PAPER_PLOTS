import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import xarray as xr
import os

file_z='maps_runoff_2200.npz'

fig, axs = plt.subplots(nrows=1,ncols=1,figsize=(18.0,18.0))

# MAR masks:
grd=xr.open_dataset('RCM_ice_regrid_04000m.nc2',decode_cf=False)
msk_ice = grd.AIS.where( (grd.AIS>0.) )* 0 + 1

# Basin masks:
basin=xr.open_dataset('/data/njourdain/DATA_ISMIP6/Mask_Iceshelf_IMBIE_4km.nc')

# Topo (just for the plot):
topo=xr.open_dataset('/data/njourdain/DATA_ISMIP6/bedmap2_8km.nc')
topo2=xr.open_dataset('/data/njourdain/DATA_ISMIP6/BedMachineAntarctica_2020-07-15_v02_8km.nc')

# Previously calculated years of hydrofracturing potential:
# First for the 16-model weighted ensemble until 2100 :
zyears=np.load('runoff_dates_2100_new.npz')
threshold = zyears['threshold']
timewin = zyears['timewin']
print('threshold = ',threshold,'kg/m2/yr ;   timewin = ',timewin,' years')
year_bas_ssp126_pct05 = np.array(zyears['year_bas_ssp126_pct05'])
year_bas_ssp126_pct17 = np.array(zyears['year_bas_ssp126_pct17'])
year_bas_ssp126_pct83 = np.array(zyears['year_bas_ssp126_pct83'])
year_bas_ssp126_pct95 = np.array(zyears['year_bas_ssp126_pct95'])
year_bas_ssp245_pct05 = np.array(zyears['year_bas_ssp245_pct05'])
year_bas_ssp245_pct17 = np.array(zyears['year_bas_ssp245_pct17'])
year_bas_ssp245_pct83 = np.array(zyears['year_bas_ssp245_pct83'])
year_bas_ssp245_pct95 = np.array(zyears['year_bas_ssp245_pct95'])
year_bas_ssp585_pct05 = np.array(zyears['year_bas_ssp585_pct05'])
year_bas_ssp585_pct17 = np.array(zyears['year_bas_ssp585_pct17'])
year_bas_ssp585_pct83 = np.array(zyears['year_bas_ssp585_pct83'])
year_bas_ssp585_pct95 = np.array(zyears['year_bas_ssp585_pct95'])
# Second for the 8-model ensemble until 2200 :
#zyearsb=np.load('runoff_dates_2200_new.npz')
zyearsb=np.load('runoff_dates_2200_weight.npz')
thresholdb = zyearsb['threshold']
timewinb = zyearsb['timewin']
if ( ( not threshold[0] == thresholdb[0] ) | ( not timewin == timewinb ) ):
  print(' *** !!! WARNIG !!! *** CHECK CONSISTENCY OF threshold AND timewin VALUES !!!!!!! ><><><><><><><><><<><><><><><<<>>>><<<<>>>>>')
year_bas_ssp126_pct05b = np.array(zyearsb['year_bas_ssp126_pct05'])
year_bas_ssp126_pct17b = np.array(zyearsb['year_bas_ssp126_pct17'])
year_bas_ssp126_pct83b = np.array(zyearsb['year_bas_ssp126_pct83'])
year_bas_ssp126_pct95b = np.array(zyearsb['year_bas_ssp126_pct95'])
year_bas_ssp585_pct05b = np.array(zyearsb['year_bas_ssp585_pct05'])
year_bas_ssp585_pct17b = np.array(zyearsb['year_bas_ssp585_pct17'])
year_bas_ssp585_pct83b = np.array(zyearsb['year_bas_ssp585_pct83'])
year_bas_ssp585_pct95b = np.array(zyearsb['year_bas_ssp585_pct95'])


model  = [ 'ACCESS-CM2', 'ACCESS-ESM1-5', 'CanESM5' , 'CESM2-WACCM', 'GISS-E2-1-H', 'IPSL-CM6A-LR',  'MRI-ESM2-0', 'UKESM1-0-LL' ]
member = [ 'r1i1p1f1'  , 'r1i1p1f1'     , 'r1i1p1f1', 'r1i1p1f1'   , 'r1i1p1f2'   , 'r1i1p1f1'    ,  'r1i1p1f1'  , 'r4i1p1f2'    ]
Nmod = np.size(model)

map_2181_2200_585 = np.zeros(np.shape(msk_ice))
nn585 = 0

############################################################

if not os.path.exists(file_z): 


   for kmod in np.arange(Nmod):
      
     print(model[kmod])
        
     file_his = 'MAR-'+model[kmod]+'_aru_1980-2014_histo_regrid_04000m.nc'
     if not os.path.exists(file_his):
       file_his = 'MAR-'+model[kmod]+'-'+member[kmod]+'_aru_1980-2014_histo_regrid_04000m_FROM_6_MODELS.nc'
     if ( model[kmod] == 'UKESM1-0-LL' ):
       file_his = 'MAR-'+model[kmod]+'-'+member[kmod]+'_aru_1980-2014_histo_regrid_04000m_FROM_UKESM1-0-LL-r1i1p1f2-histo.nc'
      
     file_ssp126 = 'MAR-'+model[kmod]+'_aru_2015-2200_ssp126_regrid_04000m.nc'
     if not os.path.exists(file_ssp126):
       file_ssp126 = 'MAR-'+model[kmod]+'_aru_2015-2200_ssp126_regrid_04000m_FROM_ssp585.nc'
       if not os.path.exists(file_ssp126):
         file_ssp126 = 'MAR-'+model[kmod]+'-'+member[kmod]+'_aru_2015-2200_ssp126_regrid_04000m_MERGED.nc'
        
     file_ssp585 = 'MAR-'+model[kmod]+'_aru_2015-2200_ssp585_regrid_04000m.nc'
     if not os.path.exists(file_ssp585):
       file_ssp585 = 'MAR-'+model[kmod]+'-'+member[kmod]+'_aru_2015-2200_ssp585_regrid_04000m_MERGED.nc'
        
     dd585=xr.open_dataset(file_ssp585,decode_cf=False)
     map_2181_2200_585 = map_2181_2200_585 + dd.aru.isel(time=slice(165,186)).mean(dim=["time"]).values
     nn585 = nn585 + 1
  
   map_2181_2200_585 = map_2181_2200_585 / nn585

   np.savez(file_z,map_2181_2200_585 = map_2181_2200_585)

else:

   zz=np.load(file_z)
   map_2181_2200 = zz['map_2181_2200']

map_2181_2200 = map_2181_2200 * msk_ice.values * 365.25 * 86400      # mm w. eq. / year

##########################################################################

# barycenters of individual ice shelves and corresponding ocean basin:
Nbasin = basin.ID.size
xbisf = np.zeros((Nbasin)) * np.nan
ybisf = np.zeros((Nbasin)) * np.nan
xboce = np.zeros((Nbasin)) * np.nan
yboce = np.zeros((Nbasin)) * np.nan
xbint = np.zeros((Nbasin)) * np.nan # intermediate point for 2-piece segments
ybint = np.zeros((Nbasin)) * np.nan
x2d, y2d = np.meshgrid(basin.x.values, basin.y.values, indexing='xy')
for kbasin in np.arange(Nbasin):
  mskisf = ( grd.ICE_MAR - grd.GROUND ) * ( basin.Iceshelf_extrap.where( (basin.Iceshelf_extrap == basin.ID.values[kbasin]) ) * 0 + 1) * grd.af2
  mskgrd = grd.GROUND * ( basin.Iceshelf_extrap.where( (basin.Iceshelf_extrap == basin.ID.values[kbasin]) ) * 0 + 1) * grd.af2
  mskoce = ( 1. - grd.AIS ) * ( basin.Iceshelf_extrap.where( (basin.Iceshelf_extrap == basin.ID.values[kbasin]) ) * 0 + 1) * grd.af2
  tmpisf = np.nansum(mskisf.values)
  tmpgrd = np.nansum(mskgrd.values)
  tmpoce = np.nansum(mskoce.values)
  # Limit of 1500 km2 for the ice shelves and 15000 km2 for the grounded part + a few exceptions
  if ( (  ( (tmpisf*4.0**2 > 1500 ) & (tmpgrd*4.0**2 > 15000 ) ) \
        | ( basin.NAME.values[kbasin] == 'LarsenA' )   \
        | ( basin.NAME.values[kbasin] == 'LarsenB' )   \
        | ( basin.NAME.values[kbasin] == 'Wordie'  )   \
        | ( basin.NAME.values[kbasin] == 'Venable' )   \
        | ( basin.NAME.values[kbasin] == 'Cosgrove') ) \
      & ( not basin.NAME.values[kbasin] == '' ) ):
    xbisf[kbasin] = np.nansum(mskisf.values * x2d) / tmpisf
    ybisf[kbasin] = np.nansum(mskisf.values * y2d) / tmpisf
    xboce[kbasin] = np.nansum(mskoce.values * x2d) / tmpoce
    yboce[kbasin] = np.nansum(mskoce.values * y2d) / tmpoce
    xbint[kbasin] = xbisf[kbasin]
    ybint[kbasin] = ybisf[kbasin]
    # Manual corrections:
    if ( basin.NAME.values[kbasin] == 'Wordie' ):
      xbint[kbasin] = xbisf[kbasin] - 0.08e6
      xboce[kbasin] = -2.5e6
      yboce[kbasin] =  1.1e6
    elif ( basin.NAME.values[kbasin] == 'Wilkins' ):
      yboce[kbasin] = yboce[kbasin] + 0.25e6
    elif ( basin.NAME.values[kbasin] == 'George_VI' ):
      xbisf[kbasin] = -2.e6
      ybisf[kbasin] =  0.77e6
      xboce[kbasin] = -2.4e6
      yboce[kbasin] =  0.80e6
      xbint[kbasin] = xbisf[kbasin]
      ybint[kbasin] = ybisf[kbasin]
    elif ( basin.NAME.values[kbasin] == 'Venable' ):  
      xboce[kbasin] = -2.5e6 
      yboce[kbasin] = ybisf[kbasin]+0.05e6
    elif ( basin.NAME.values[kbasin] == 'Abbot' ): 
      yboce[kbasin] = ybisf[kbasin]+0.05e6
    elif ( basin.NAME.values[kbasin] == 'Cosgrove' ): 
      xbint[kbasin] = xbisf[kbasin]
      ybint[kbasin] = ybisf[kbasin] - 0.12e6
      xboce[kbasin] = -2.52e6
      yboce[kbasin] = ybint[kbasin] - 0.05e6
    elif ( basin.NAME.values[kbasin] == 'Pine_Island' ):
      xbint[kbasin] = xbisf[kbasin] - 0.05e6
      ybint[kbasin] = ybisf[kbasin] - 0.15e6
      xboce[kbasin] = xbisf[kbasin] - 0.33e6
      yboce[kbasin] = ybint[kbasin] - 0.08e6
    elif ( basin.NAME.values[kbasin] == 'Thwaites' ):
      xbisf[kbasin] = xbisf[kbasin] - 0.02e6
      xbint[kbasin] = -2.0e6
      ybint[kbasin] = -0.9e6
      xboce[kbasin] = -2.5e6
      yboce[kbasin] = ybint[kbasin] - 0.17e6
    elif ( basin.NAME.values[kbasin] == 'Dotson/Philbin_Inlet' ):
      xboce[kbasin] = xboce[kbasin] - 0.25e6
      yboce[kbasin] = yboce[kbasin] - 1.0e6
    elif ( basin.NAME.values[kbasin] == 'Crosson' ):
      xbint[kbasin] = xboce[kbasin] - 0.02e6
      ybint[kbasin] = yboce[kbasin] + 0.02e6
      xboce[kbasin] = -2.1e6
      yboce[kbasin] = -1.4e6
    elif ( basin.NAME.values[kbasin] == 'Getz' ):
      xboce[kbasin] = xbisf[kbasin] - 0.02e6
      yboce[kbasin] = yboce[kbasin] - 0.75 * (yboce[kbasin] - ybisf[kbasin])
    elif ( basin.NAME.values[kbasin] == 'Sulzberger' ):
      xboce[kbasin] = xboce[kbasin] - 2 * (xboce[kbasin] - xbisf[kbasin])
    elif ( basin.NAME.values[kbasin] == 'Ross_West' ):
      xboce[kbasin] = -0.6e6
      yboce[kbasin] = -0.5e6
    elif ( basin.NAME.values[kbasin] == 'Ross_East' ):
      xboce[kbasin] = 0.5e6
      yboce[kbasin] = -0.5e6
    elif ( basin.NAME.values[kbasin] == 'Nordenskjold/Marin/HarbordGlacier/Cheetham/GeikieInlet' ):
      xbint[kbasin] = 0.
      xboce[kbasin] = -0.45e6
      yboce[kbasin] = -1.8e6
    elif ( basin.NAME.values[kbasin] == 'Drygalski' ):
      xbint[kbasin] = xboce[kbasin]
      ybint[kbasin] = yboce[kbasin]
      xboce[kbasin] =  0.1e6
      yboce[kbasin] = -1.7e6      
    elif ( basin.NAME.values[kbasin] == 'Nansen' ):
      xbint[kbasin] = xboce[kbasin] - 0.15e6
      ybint[kbasin] = ybisf[kbasin] - 0.05e6
      xboce[kbasin] =  0.05e6
      yboce[kbasin] = -2.1e6
    elif ( basin.NAME.values[kbasin] == 'Rennick' ):
      xboce[kbasin] = xbisf[kbasin] - 0.05e6
    elif ( basin.NAME.values[kbasin] == 'Cook' ):
      xboce[kbasin] = xbisf[kbasin] - 0.05e6
      yboce[kbasin] = yboce[kbasin] - 0.13e6
    elif ( basin.NAME.values[kbasin] == 'Ninnis' ):
      xboce[kbasin] = xbisf[kbasin] - 0.05e6
      yboce[kbasin] = -2.8e6
    elif ( basin.NAME.values[kbasin] == 'Mertz' ):
      xboce[kbasin] = xbisf[kbasin] + 0.02e6
    elif ( basin.NAME.values[kbasin] == 'WattBay/Zelee/Astrolabe/Liotard/Francais/Marret/Commandant_Charcot//PourquoiPas' ):
      yboce[kbasin] = -2.9e6
    elif ( basin.NAME.values[kbasin] == 'Dibble' ):
      xboce[kbasin] = 2.4e6
    elif ( basin.NAME.values[kbasin] == 'May_Glacier/Morse/Sandford' ):
      yboce[kbasin] = yboce[kbasin] -0.08e6
    elif ( basin.NAME.values[kbasin] == 'Moscow_University' ):
      xboce[kbasin] = 2.7e6
    elif ( basin.NAME.values[kbasin] == 'Totten' ):
      xboce[kbasin] = 2.8e6
      yboce[kbasin] = ybisf[kbasin] - 0.05e6
    elif ( basin.NAME.values[kbasin] == 'Shackleton' ):
      yboce[kbasin] = ybisf[kbasin] - ( yboce[kbasin] - ybisf[kbasin] )
    elif ( basin.NAME.values[kbasin] == 'Helen' ):
      yboce[kbasin] = ybisf[kbasin] + 0.01e6
    elif ( basin.NAME.values[kbasin] == 'West' ):
      xboce[kbasin] = 3.2e6
      yboce[kbasin] = ybisf[kbasin] - 0.05e6
    elif ( basin.NAME.values[kbasin] == 'Publications' ):
      xboce[kbasin] = 2.9e6 
      yboce[kbasin] = 1.0e6 
    elif ( basin.NAME.values[kbasin] == 'Amery' ):
      xbint[kbasin] = xboce[kbasin]
      ybint[kbasin] = ybisf[kbasin] + 0.03e6
      xboce[kbasin] = 2.8e6
      yboce[kbasin] = 1.4e6
    elif ( basin.NAME.values[kbasin] == 'Utsikkar/Mulebreen/Cirque_Fjord/Hoseason/Rund_Bay' ):
      xboce[kbasin] = 2.4e6
      yboce[kbasin] = 1.8e6
    elif ( basin.NAME.values[kbasin] == 'Zubchatyy/Porter/Myers' ):
      yboce[kbasin] = yboce[kbasin] +0.2e6
    elif ( basin.NAME.values[kbasin] == 'Hannan/Telen/Skallen' ):
      xbint[kbasin] = xboce[kbasin] + 0.1e6
      ybint[kbasin] = ybint[kbasin] + 0.1e6
      xboce[kbasin] = 2.6e6
      yboce[kbasin] = 2.7e6
    elif ( basin.NAME.values[kbasin] == 'Prince_Harald' ):
      yboce[kbasin] = yboce[kbasin] + 0.2e6
    elif ( basin.NAME.values[kbasin] == 'Baudouin' ):
      xboce[kbasin] = 1.9e6
      yboce[kbasin] = 2.7e6
    elif ( basin.NAME.values[kbasin] == 'Borchgrevink' ):
      yboce[kbasin] = yboce[kbasin] - 0.1e6
    elif ( basin.NAME.values[kbasin] == 'Lazarev' ):
      xboce[kbasin] = xboce[kbasin] + 0.4e6
    elif ( basin.NAME.values[kbasin] == 'Nivl' ):
      ybint[kbasin] = 2.6e6
      xboce[kbasin] = xbisf[kbasin] + 0.05e6
      yboce[kbasin] = 2.7e6
    elif ( basin.NAME.values[kbasin] == 'Vigrid' ):
      xbint[kbasin] = xbisf[kbasin] + 0.05e6
      ybint[kbasin] = 2.5e6
      xboce[kbasin] = xbint[kbasin] - 0.05e6
      yboce[kbasin] = 2.7e6
    elif ( basin.NAME.values[kbasin] == 'Fimbul' ):
      xbisf[kbasin] = xbisf[kbasin] - 0.05e6
      xbint[kbasin] = xbint[kbasin] - 0.05e6
      xboce[kbasin] = xboce[kbasin] - 0.15e6
      ybisf[kbasin] = ybisf[kbasin]
      ybint[kbasin] = ybint[kbasin]
      yboce[kbasin] = 2.7e6
    elif ( basin.NAME.values[kbasin] == 'Jelbart' ):
      yboce[kbasin] = yboce[kbasin] + 0.1e6
    elif ( basin.NAME.values[kbasin] == 'Ekstrom' ):
      xbint[kbasin] = xbisf[kbasin] - 0.35e6
      ybint[kbasin] = ybisf[kbasin] + 0.04e6
      xboce[kbasin] = xboce[kbasin] - 0.20e6
      yboce[kbasin] = 2.7e6
    elif ( basin.NAME.values[kbasin] == 'Riiser-Larsen' ):
      xboce[kbasin] = xboce[kbasin] + 0.30e6
      yboce[kbasin] = yboce[kbasin] + 0.05e6
    elif ( basin.NAME.values[kbasin] == 'Brunt_Stancomb' ):
      xboce[kbasin] = xboce[kbasin] + 0.30e6
      yboce[kbasin] = yboce[kbasin] + 0.20e6
    elif ( basin.NAME.values[kbasin] == 'Dawson_Lambton/Hayes_Coats_Coast' ):
      xboce[kbasin] = xboce[kbasin] + 0.08e6
      yboce[kbasin] = yboce[kbasin] + 0.15e6
    elif ( basin.NAME.values[kbasin] == 'Filchner' ):
      xboce[kbasin] =  0.e0
      yboce[kbasin] =  1.0e6
    elif ( basin.NAME.values[kbasin] == 'Ronne' ):
      xbisf[kbasin] = xbisf[kbasin] + 0.1e6
      xbint[kbasin] = xbint[kbasin] + 0.1e6
      xboce[kbasin] = xboce[kbasin] + 0.1e6
      ybisf[kbasin] = ybisf[kbasin] + 0.1e6
      ybint[kbasin] = ybint[kbasin] + 0.1e6
    elif ( basin.NAME.values[kbasin] == 'LarsenD' ):
      xboce[kbasin] = xboce[kbasin] - 0.1e6
      yboce[kbasin] = yboce[kbasin] + 0.1e6
      xbint[kbasin] = xboce[kbasin] - 0.05e6
      ybint[kbasin] = yboce[kbasin] - 0.05e6
    elif ( basin.NAME.values[kbasin] == 'LarsenC' ):
      xboce[kbasin] = xboce[kbasin] - 0.05e6
      yboce[kbasin] = 1.8e6
    elif ( basin.NAME.values[kbasin] == 'LarsenB' ):
      xboce[kbasin] = xboce[kbasin] + 0.2e6
      yboce[kbasin] = 2.2e6
      xbint[kbasin] = xboce[kbasin] - 0.1e6
      ybint[kbasin] = yboce[kbasin] - 0.1e6
    elif ( basin.NAME.values[kbasin] == 'LarsenA' ):
      xbint[kbasin] = xbisf[kbasin]
      ybint[kbasin] = 2.60e6
      xboce[kbasin] = xboce[kbasin] + 0.2e6
      yboce[kbasin] = 2.70e6
    # Correction of names:
    if ( basin.NAME.values[kbasin] == 'Utsikkar/Mulebreen/Cirque_Fjord/Hoseason/Rund_Bay' ):
      basin.NAME.values[kbasin] = 'Utsikkar, Mulebreen, Cirque Fjord,\nHoseason, Rund Bay'
    elif ( basin.NAME.values[kbasin] == 'Nordenskjold/Marin/HarbordGlacier/Cheetham/GeikieInlet' ):
      basin.NAME.values[kbasin] = 'Nordenskjöld, Marin, Harbord,\nCheetham, Geikie Inlet'
    elif ( basin.NAME.values[kbasin] == 'WattBay/Zelee/Astrolabe/Liotard/Francais/Marret/Commandant_Charcot//PourquoiPas' ):
      basin.NAME.values[kbasin] = 'Watt Bay, Zélée, Astrolabe, Liotard, Barré,\nFrançais, Marret, Comdt Charcot, Pourquoi Pas'
    elif ( basin.NAME.values[kbasin] == 'Baudouin' ):
      basin.NAME.values[kbasin] = 'Roi Baudouin'
    elif ( basin.NAME.values[kbasin] == 'Dotson/Philbin_Inlet' ):
      basin.NAME.values[kbasin] = 'Dotson'
    elif ( basin.NAME.values[kbasin] == 'Ekstrom' ):
      basin.NAME.values[kbasin] = 'Ekström'
    elif ( basin.NAME.values[kbasin] == 'Dawson_Lambton/Hayes_Coats_Coast' ):
      basin.NAME.values[kbasin] = 'Dawson-\nLambton,\nHayes'
    basin.NAME.values[kbasin] = basin.NAME.values[kbasin].replace("_"," ")
    basin.NAME.values[kbasin] = basin.NAME.values[kbasin].replace("/",", ")
    #-
    print(basin.ID.values[kbasin], basin.NAME.values[kbasin], tmpisf*4.0**2, tmpoce*4.0**2, xboce[kbasin], yboce[kbasin])

##########################################################################
# PLOT :

#----------
# Colorbar:
cbar_range = np.arange(0.,4200.,200.)

#----------
# Defining colormap:

# moving the zero of colorbar
# NB: modify the Ncool to Nwarm ratio (total=256) to place zero as desired.
Ncool=int(256*(-np.amin(cbar_range)/(np.amax(cbar_range)-np.amin(cbar_range))))
Nwarm=256-Ncool
col = cm.get_cmap('PuOr', 256)
#tmp1 = col(np.linspace(0.47, 1.00, Ncool)) # decrease first number to have more white in the middle light-blue colors
#tmp2 = col(np.linspace(0.00, 0.51, Nwarm)) # increase second number to have more white in the middle light-yellow colors
#newcolors = np.append(tmp1[::-1,:],tmp2[::-1,:],axis=0)
newcolors = col(np.linspace(0.47, 1.00,256))
cmap_new = ListedColormap(newcolors)

cax=fig.add_axes([0.52, 0.02, 0.43, 0.015]) # color bar
im0=axs.contourf(grd.x,grd.y,map_2181_2200,cbar_range,cmap=cmap_new,extend='max')
cbar0=fig.colorbar(im0, cax=cax, orientation="horizontal")
cbar0.ax.tick_params(labelsize=14)
cbar0.set_label(r'Anomaly of liquid water production in excess (kg m$^{-2}$ yr$^{-1}$'+')\n'+ 'in 2181-2200 w.r.t. 1995-2014 (SSP5-8.5)',labelpad=-80,fontsize=14)

axs.contour(topo.x,topo.y,topo.surface.values,[1000,2000,3000,4000],colors='darkgrey',linewidths=1)
axs.contour(topo2.x,topo2.y,topo2.mask.values,[0.5],colors='black',linewidths=1)
axs.contour(basin.x,basin.y,basin.Iceshelf,np.arange(0.5,200.5,1),colors='black',linewidths=2.5)

for kbasin in np.arange(Nbasin):

  if ( not np.isnan(xboce[kbasin]) ):

    #-----------------------------------------
    # combining the 2 ensembles:
    # 8-member ensemble if year>2100 in 16-member weighted ensemble
    for scenar in ['ssp126', 'ssp585']:
       for pct in ['05', '17', '83', '95']:
          if ( eval("year_bas_"+scenar+"_pct"+pct+"[kbasin]") == 2500 ):
             exec("year_bas_"+scenar+"_pct"+pct+"[kbasin] = np.max([2101, year_bas_"+scenar+"_pct"+pct+"b[kbasin]])")
 
    #-----------------------------------------
    # Ice Shelf Labels :

    dx = (xboce[kbasin]-xbint[kbasin]) / np.sqrt( (xboce[kbasin]-xbint[kbasin])**2 + (yboce[kbasin]-ybint[kbasin])**2 ) * 1.e4
    dy = (yboce[kbasin]-ybint[kbasin]) / np.sqrt( (xboce[kbasin]-xbint[kbasin])**2 + (yboce[kbasin]-ybint[kbasin])**2 ) * 1.e4

    if ( dx < 0. ):
      hal = 'right'
    else:
      hal = 'left'
    #--
    if ( dy < 0. ):
      val = 'top'
    else:
      val = 'bottom'
    #--
    if (  ( basin.NAME[kbasin].values == 'Nordenskjöld, Marin, Harbord,\nCheetham, Geikie Inlet' ) \
        | ( basin.NAME[kbasin].values == 'Watt Bay, Zélée, Astrolabe, Liotard, Barré,\nFrançais, Marret, Comdt Charcot, Pourquoi Pas' ) ):
      ddyy = 8.0e4
    else:
      ddyy = 0.e0
    #--
    if ( (year_bas_ssp126_pct83[kbasin] == 1500) & (year_bas_ssp126_pct17[kbasin] == 1500) ):
      textssp126='always'
    elif ( (year_bas_ssp126_pct83[kbasin] == 2500) & (year_bas_ssp126_pct17[kbasin] == 2500) ):
      textssp126='never'
    elif (year_bas_ssp126_pct83[kbasin] == 1500):
      textssp126='alw.-'+year_bas_ssp126_pct17[kbasin].astype('str')
    elif (year_bas_ssp126_pct17[kbasin] == 2500):
      textssp126=year_bas_ssp126_pct83[kbasin].astype('str')+'-nev.'
    else:
      textssp126=year_bas_ssp126_pct83[kbasin].astype('str')+'-'+year_bas_ssp126_pct17[kbasin].astype('str')
    #--
    if ( (year_bas_ssp585_pct83[kbasin] == 1500) & (year_bas_ssp585_pct17[kbasin] == 1500) ):
      textssp585='always'
    elif ( (year_bas_ssp585_pct83[kbasin] == 2500) & (year_bas_ssp585_pct17[kbasin] == 2500) ):
      textssp585='never'
    elif (year_bas_ssp585_pct83[kbasin] == 1500):
      textssp585='alw.-'+year_bas_ssp585_pct17[kbasin].astype('str')
    elif (year_bas_ssp585_pct17[kbasin] == 2500):
      textssp585=year_bas_ssp585_pct83[kbasin].astype('str')+'-nev.'
    else:
      textssp585=year_bas_ssp585_pct83[kbasin].astype('str')+'-'+year_bas_ssp585_pct17[kbasin].astype('str')
    #--
    if (year_bas_ssp126_pct17[kbasin] <= 2014):
       coltex = 'darkred' # Full likely range before 2015
    elif ( (year_bas_ssp126_pct83[kbasin] <= 2014) | (year_bas_ssp585_pct83[kbasin] <= 2014) ):
       coltex = 'red'     # Likely range starting before 2015
    elif (year_bas_ssp585_pct83[kbasin] <= 2050):
       coltex = 'orchid'  # likely range starting before 2050 in at least one SSP
    elif (year_bas_ssp585_pct83[kbasin] <= 2100):
       coltex = 'blueviolet'  # likely range starting before 2100 under ssp585
    else:
       coltex = 'darkblue' # likely range starting after 2100 under ssp585
    #--
    axs.plot([xbisf[kbasin],xbint[kbasin],xboce[kbasin]],[ybisf[kbasin],ybint[kbasin],yboce[kbasin]],color=coltex)
    axs.text(xboce[kbasin]+dx,yboce[kbasin]+dy,basin.NAME.values[kbasin],fontsize=12,color=coltex,ha=hal,va=val,fontweight='bold')
    axs.text(xboce[kbasin]+dx,yboce[kbasin]+dy-8.0e4-ddyy,textssp126,fontsize=12,color=coltex,ha=hal,va=val,fontstyle='italic')
    axs.text(xboce[kbasin]+dx,yboce[kbasin]+dy-16.e4-ddyy,textssp585,fontsize=12,color=coltex,ha=hal,va=val)

# Legend for color code:
axs.text(-2.80e6,-2.3e6,'Color code :',color='black',fontsize=13)
axs.text(-2.80e6,-2.4e6,'Full likely range before 2015',color='darkred',fontsize=13,fontweight='bold')
axs.text(-2.80e6,-2.5e6,'Likely range starting before 2015',color='red',fontsize=13,fontweight='bold')
axs.text(-2.80e6,-2.6e6,'Likely range starting before 2050 in at least one SSP',color='orchid',fontsize=13,fontweight='bold')
axs.text(-2.80e6,-2.7e6,'Likely range starting before 2100 in SSP5-8.5',color='blueviolet',fontsize=13,fontweight='bold')
axs.text(-2.80e6,-2.8e6,'Likely range starting after 2100 in SSP5-8.5',color='darkblue',fontsize=13,fontweight='bold')

# Additional panel on the evolution of number of ice shelves above the threshold :
yyyy=np.arange(1850,2201,1)
Ny=np.size(yyyy)
nb_isf_ssp126_pct83 = np.zeros((Ny))
nb_isf_ssp126_pct17 = np.zeros((Ny))
nb_isf_ssp585_pct83 = np.zeros((Ny))
nb_isf_ssp585_pct17 = np.zeros((Ny))
for ky in np.arange(Ny):
   nb_isf_ssp126_pct83[ky] = np.sum( year_bas_ssp126_pct83 <= yyyy[ky] )
   nb_isf_ssp126_pct17[ky] = np.sum( year_bas_ssp126_pct17 <= yyyy[ky] )
   nb_isf_ssp585_pct83[ky] = np.sum( year_bas_ssp585_pct83 <= yyyy[ky] )
   nb_isf_ssp585_pct17[ky] = np.sum( year_bas_ssp585_pct17 <= yyyy[ky] )
delta = nb_isf_ssp126_pct83[0]
nb_isf_ssp126_pct83 = nb_isf_ssp126_pct83 - delta + 6
nb_isf_ssp126_pct17 = nb_isf_ssp126_pct17 - delta + 6
nb_isf_ssp585_pct83 = nb_isf_ssp585_pct83 - delta + 6
nb_isf_ssp585_pct17 = nb_isf_ssp585_pct17 - delta + 6
llax=fig.add_axes([0.17, 0.02, 0.32, 0.10])
llax.fill_between(yyyy,nb_isf_ssp126_pct83,nb_isf_ssp126_pct17,color='cornflowerblue',alpha=0.2)
llax.fill_between(yyyy,nb_isf_ssp585_pct83,nb_isf_ssp585_pct17,color='darkred',alpha=0.2)
llax.text(1852,50,'Number of ice shelves above threshold',fontsize=14)
llax.text(2115,40,'SSP5-8.5',color='darkred',fontsize=15)
llax.text(2140,15,'SSP1-2.6',color='cornflowerblue',fontsize=15)
llax.tick_params(axis='both', labelsize=14)
llax.set_xlim([1850,2200])
llax.set_ylim([0,np.max(nb_isf_ssp585_pct17)])

xc=np.mean(axs.get_xlim())
yu=axs.get_ylim()[1]
axs.text(xc,0.98*yu,'Likely emergence of surface conditions necessary for hydrofracturing',fontsize=20,fontweight='bold',ha='center')
axs.set_axis_off()

##########################################################################

fig.savefig("map_runoff_projection_new.pdf")
