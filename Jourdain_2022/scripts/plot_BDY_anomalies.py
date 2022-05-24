import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

file_TS_p = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MAR/bdyT_tra_climato.nc'
file_TS_f = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MARrcp85/bdyT_tra_climato.nc'
file_Ubtp_p = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MAR/bdyU_u2d_climato.nc'
file_Ubtp_f = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MARrcp85/bdyU_u2d_climato.nc'
file_Ubcl_p = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MAR/bdyU_u3d_climato.nc'
file_Ubcl_f = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MARrcp85/bdyU_u3d_climato.nc'
file_Vbtp_p = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MAR/bdyV_u2d_climato.nc'
file_Vbtp_f = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MARrcp85/bdyV_u2d_climato.nc'
file_Vbcl_p = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MAR/bdyV_u3d_climato.nc'
file_Vbcl_f = '../OUTPUT_AMU/nemo_AMUXL12_GNJ002_BM02MARrcp85/bdyV_u3d_climato.nc'

nc_TS_p = xr.open_dataset(file_TS_p,decode_cf=False).squeeze('yb')
nc_TS_f = xr.open_dataset(file_TS_f,decode_cf=False).squeeze('yb')
nc_Ubtp_p = xr.open_dataset(file_Ubtp_p,decode_cf=False).squeeze('yb')
nc_Ubtp_f = xr.open_dataset(file_Ubtp_f,decode_cf=False).squeeze('yb')
nc_Ubcl_p = xr.open_dataset(file_Ubcl_p,decode_cf=False).squeeze('yb')
nc_Ubcl_f = xr.open_dataset(file_Ubcl_f,decode_cf=False).squeeze('yb')
nc_Vbtp_p = xr.open_dataset(file_Vbtp_p,decode_cf=False).squeeze('yb')
nc_Vbtp_f = xr.open_dataset(file_Vbtp_f,decode_cf=False).squeeze('yb')
nc_Vbcl_p = xr.open_dataset(file_Vbcl_p,decode_cf=False).squeeze('yb')
nc_Vbcl_f = xr.open_dataset(file_Vbcl_f,decode_cf=False).squeeze('yb')

Stmp = nc_TS_p.vosaline.isel(time_counter=0)
msk = Stmp.where( (Stmp == 0.), np.nan)

latT = nc_TS_p.nav_lat
lonT = nc_TS_p.nav_lon
latU = nc_Ubtp_p.nav_lat
lonU = nc_Ubtp_p.nav_lon
latV = nc_Vbtp_p.nav_lat
lonV = nc_Vbtp_p.nav_lon

# future minus present anomalies:
T_anom = nc_TS_f.votemper.mean(axis=0) - nc_TS_p.votemper.mean(axis=0)
S_anom = nc_TS_f.vosaline.mean(axis=0) - nc_TS_p.vosaline.mean(axis=0)
# anomalies combining barotropic and baroclinic components:
tmpU = nc_Ubtp_f.vobtcrtx.mean(axis=0) - nc_Ubtp_p.vobtcrtx.mean(axis=0)
tmpV = nc_Vbtp_f.vobtcrty.mean(axis=0) - nc_Vbtp_p.vobtcrty.mean(axis=0)
U_anom = xr.DataArray( np.repeat(tmpU.values[ np.newaxis, :],75, axis=0), dims=['depthu','xbU'] )
V_anom = xr.DataArray( np.repeat(tmpV.values[ np.newaxis, :],75, axis=0), dims=['depthv','xbV'] )
U_anom = U_anom + nc_Ubcl_f.vozocrtx.mean(axis=0) - nc_Ubcl_p.vozocrtx.mean(axis=0)
V_anom = V_anom + nc_Vbcl_f.vomecrty.mean(axis=0) - nc_Vbcl_p.vomecrty.mean(axis=0)

#----------------------------------------------------------------------
# reordering indices:

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

k1ET=0    ; k2ET=564  ; nET=k2ET-k1ET+1
k1WT=565  ; k2WT=1128 ; nWT=k2WT-k1WT+1
k1NT=1129 ; k2NT=1812 ; nNT=k2NT-k1NT+1

k1EU=0    ; k2EU=564  ; nEU=k2EU-k1EU+1
k1WU=565  ; k2WU=1128 ; nWU=k2WU-k1WU+1
k1NU=1129 ; k2NU=1811 ; nNU=k2NU-k1NU+1

k1EV=0    ; k2EV=565  ; nEV=k2EV-k1EV+1
k1WV=566  ; k2WV=1130 ; nWV=k2WV-k1WV+1
k1NV=1131 ; k2NV=1814 ; nNV=k2NV-k1NV+1

msk_reshaped = np.zeros(msk.shape)
msk_reshaped[:,0:nWT] = msk.isel(xbT=slice(k1WT,k2WT+1)).values
msk_reshaped[:,nWT:nWT+nNT] = msk.isel(xbT=slice(k1NT,k2NT+1)).values
msk_reshaped[:,nWT+nNT:nWT+nNT+nET] = np.flip( msk.isel(xbT=slice(k1ET,k2ET+1)).values, axis=1 )

T_anom_reshaped = np.zeros(T_anom.shape)
T_anom_reshaped[:,0:nWT] = T_anom.isel(xbT=slice(k1WT,k2WT+1)).values
T_anom_reshaped[:,nWT:nWT+nNT] = T_anom.isel(xbT=slice(k1NT,k2NT+1)).values
T_anom_reshaped[:,nWT+nNT:nWT+nNT+nET] = np.flip( T_anom.isel(xbT=slice(k1ET,k2ET+1)).values, axis=1 )
#
S_anom_reshaped = np.zeros(S_anom.shape)
S_anom_reshaped[:,0:nWT] = S_anom.isel(xbT=slice(k1WT,k2WT+1)).values
S_anom_reshaped[:,nWT:nWT+nNT] = S_anom.isel(xbT=slice(k1NT,k2NT+1)).values
S_anom_reshaped[:,nWT+nNT:nWT+nNT+nET] = np.flip( S_anom.isel(xbT=slice(k1ET,k2ET+1)).values, axis=1 )
#
latT_reshaped = np.zeros(latT.shape)
latT_reshaped[0:nWT] = latT.isel(xbT=slice(k1WT,k2WT+1)).values
latT_reshaped[nWT:nWT+nNT] = latT.isel(xbT=slice(k1NT,k2NT+1)).values
latT_reshaped[nWT+nNT:nWT+nNT+nET] = np.flip( latT.isel(xbT=slice(k1ET,k2ET+1)).values )
#
lonT_reshaped = np.zeros(lonT.shape)
lonT_reshaped[0:nWT] = lonT.isel(xbT=slice(k1WT,k2WT+1)).values
lonT_reshaped[nWT:nWT+nNT] = lonT.isel(xbT=slice(k1NT,k2NT+1)).values
lonT_reshaped[nWT+nNT:nWT+nNT+nET] = np.flip( lonT.isel(xbT=slice(k1ET,k2ET+1)).values )
#
klatWT=np.zeros((4))
klatET=np.zeros((4))
kk=0
for xlat in [-75.,-70.,-65.,-60.]:
  klatWT[kk] = find_nearest(latT_reshaped[0:nWT],xlat)
  klatET[kk] = nWT+nNT+find_nearest(latT_reshaped[nWT+nNT:nWT+nNT+nET],xlat)
  kk=kk+1
klonNT=np.zeros((4))
kk=0
for xlon in [-135.,-120.,-105.,-90.]:
  klonNT[kk] = nWT+find_nearest(lonT_reshaped[nWT:nWT+nNT],xlon)
  kk=kk+1

U_anom_reshaped = np.zeros(U_anom.shape)
U_anom_reshaped[:,0:nWU] = U_anom.isel(xbU=slice(k1WU,k2WU+1)).values
U_anom_reshaped[:,nWU:nWU+nNU] = U_anom.isel(xbU=slice(k1NU,k2NU+1)).values
U_anom_reshaped[:,nWU+nNU:nWU+nNU+nEU] = np.flip( U_anom.isel(xbU=slice(k1EU,k2EU+1)).values, axis=1 )
#
latU_reshaped = np.zeros(latU.shape)
latU_reshaped[0:nWU] = latU.isel(xbU=slice(k1WU,k2WU+1)).values
latU_reshaped[nWU:nWU+nNU] = latU.isel(xbU=slice(k1NU,k2NU+1)).values
latU_reshaped[nWU+nNU:nWU+nNU+nEU] = np.flip( latU.isel(xbU=slice(k1EU,k2EU+1)).values )
#
lonU_reshaped = np.zeros(lonU.shape)
lonU_reshaped[0:nWU] = lonU.isel(xbU=slice(k1WU,k2WU+1)).values
lonU_reshaped[nWU:nWU+nNU] = lonU.isel(xbU=slice(k1NU,k2NU+1)).values
lonU_reshaped[nWU+nNU:nWU+nNU+nEU] = np.flip( lonU.isel(xbU=slice(k1EU,k2EU+1)).values )
#
klatWU=np.zeros((4))
klatEU=np.zeros((4))
kk=0
for xlat in [-75.,-70.,-65.,-60.]:
  klatWU[kk] = find_nearest(latU_reshaped[0:nWU],xlat)
  klatEU[kk] = nWU+nNU+find_nearest(latU_reshaped[nWU+nNU:nWU+nNU+nEU],xlat)
  kk=kk+1
klonNU=np.zeros((4))
kk=0
for xlon in [-135.,-120.,-105.,-90.]:
  klonNU[kk] = nWU+find_nearest(lonU_reshaped[nWU:nWU+nNU],xlon)
  kk=kk+1

V_anom_reshaped = np.zeros(V_anom.shape)
V_anom_reshaped[:,0:nWV] = V_anom.isel(xbV=slice(k1WV,k2WV+1)).values
V_anom_reshaped[:,nWV:nWV+nNV] = V_anom.isel(xbV=slice(k1NV,k2NV+1)).values
V_anom_reshaped[:,nWV+nNV:nWV+nNV+nEV] = np.flip( V_anom.isel(xbV=slice(k1EV,k2EV+1)).values, axis=1 )
#
latV_reshaped = np.zeros(latV.shape)
latV_reshaped[0:nWV] = latV.isel(xbV=slice(k1WV,k2WV+1)).values
latV_reshaped[nWV:nWV+nNV] = latV.isel(xbV=slice(k1NV,k2NV+1)).values
latV_reshaped[nWV+nNV:nWV+nNV+nEV] = np.flip( latV.isel(xbV=slice(k1EV,k2EV+1)).values )
#
lonV_reshaped = np.zeros(lonV.shape)
lonV_reshaped[0:nWV] = lonV.isel(xbV=slice(k1WV,k2WV+1)).values
lonV_reshaped[nWV:nWV+nNV] = lonV.isel(xbV=slice(k1NV,k2NV+1)).values
lonV_reshaped[nWV+nNV:nWV+nNV+nEV] = np.flip( lonV.isel(xbV=slice(k1EV,k2EV+1)).values )
#
klatWV=np.zeros((4))
klatEV=np.zeros((4))
kk=0
for xlat in [-75.,-70.,-65.,-60.]:
  klatWV[kk] = find_nearest(latV_reshaped[0:nWV],xlat)
  klatEV[kk] = nWV+nNV+find_nearest(latV_reshaped[nWV+nNV:nWV+nNV+nEV],xlat)
  kk=kk+1
klonNV=np.zeros((4))
kk=0
for xlon in [-135.,-120.,-105.,-90.]:
  klonNV[kk] = nWV+find_nearest(lonV_reshaped[nWV:nWV+nNV],xlon)
  kk=kk+1

#----------------------------------------------------------------------
# Plots:

fig, axs = plt.subplots(nrows=3,ncols=1,figsize=(18.0,18.0))
fig.tight_layout(pad=6.0)
axs = axs.ravel()

depmax=2000.
xxT=np.arange(np.size(latT_reshaped))
xxU=np.arange(np.size(latU_reshaped))
xxV=np.arange(np.size(latV_reshaped))

cax0=axs[0].contourf(xxT,nc_TS_p.deptht,T_anom_reshaped,np.arange(-2.25,2.5,0.25),cmap='RdBu_r')
axs[0].pcolormesh(xxT,nc_TS_p.deptht,msk_reshaped,vmin=-2,vmax=2,cmap='Greys',rasterized=True)
axs[0].plot([nWT, nWT],[depmax,0.],color='k',linewidth=1.9)
axs[0].plot([nWT+nNT, nWT+nNT],[depmax,0.],color='k',linewidth=1.9)
for kk in np.arange(np.size(klatWT)):
  axs[0].plot([klatWT[kk],klatWT[kk]],[depmax,0.],color='k',linestyle='--',linewidth=0.6)
  axs[0].plot([klatET[kk],klatET[kk]],[depmax,0.],color='k',linestyle='--',linewidth=0.6)
for kk in np.arange(np.size(klonNT)):
  axs[0].plot([klonNT[kk],klonNT[kk]],[depmax,0.],color='k',linestyle='--',linewidth=0.6)
axs[0].set_title('(a) Conservative temperature anomaly',fontsize=20,weight='bold')
axs[0].set_ylabel('Depth (m)',fontsize=18)
axs[0].invert_yaxis()
axs[0].set_ylim(depmax,0)
axs[0].set_xticks(np.append(np.append(klatWT,klonNT),np.flip(klatET)))
axs[0].set_xticklabels(['75°S','70°S','65°S','60°S','135°W','120°W','105°W','90°W','60°S','65°S','70°S','75°S'])
axs[0].tick_params(axis='x', labelsize=16)
axs[0].tick_params(axis='y', labelsize=16)
axs[0].text(klatWT[0]+90,1950,'WESTERN BOUNDARY',fontsize=12, weight='bold')
axs[0].text(klatET[3]+20,1950,'EASTERN BOUNDARY',fontsize=12, weight='bold')
axs[0].text(klonNT[0]+130,1950,'NORTHERN BOUNDARY',fontsize=12, weight='bold')
#colorbar:
cbar0=fig.colorbar(cax0,ax=axs[0],fraction=0.029, pad=0.02, ticks=np.arange(-2.,2.5,0.5))
cbar0.ax.set_title('degC',size=16)
cbar0.outline.set_linewidth(0.4)
cbar0.ax.tick_params(labelsize=16,which='both')

cax1=axs[1].contourf(xxT,nc_TS_p.deptht,S_anom_reshaped,np.arange(-0.65,0.70,0.05),cmap='RdBu_r')
axs[1].pcolormesh(xxT,nc_TS_p.deptht,msk_reshaped,vmin=-2,vmax=2,cmap='Greys',rasterized=True)
axs[1].plot([nWT, nWT],[depmax,0.],color='k',linewidth=1.9)
axs[1].plot([nWT+nNT, nWT+nNT],[depmax,0.],color='k',linewidth=1.9)
for kk in np.arange(np.size(klatWT)):
  axs[1].plot([klatWT[kk],klatWT[kk]],[depmax,0.],color='k',linestyle='--',linewidth=0.6)
  axs[1].plot([klatET[kk],klatET[kk]],[depmax,0.],color='k',linestyle='--',linewidth=0.6)
for kk in np.arange(np.size(klonNT)):
  axs[1].plot([klonNT[kk],klonNT[kk]],[depmax,0.],color='k',linestyle='--',linewidth=0.6)
axs[1].set_title('(b) Absolute salinity anomaly',fontsize=20,weight='bold')
axs[1].set_ylabel('Depth (m)',fontsize=18)
axs[1].invert_yaxis()
axs[1].set_ylim(depmax,0)
axs[1].set_xticks(np.append(np.append(klatWT,klonNT),np.flip(klatET)))
axs[1].set_xticklabels(['75°S','70°S','65°S','60°S','135°W','120°W','105°W','90°W','60°S','65°S','70°S','75°S'])
axs[1].tick_params(axis='x', labelsize=16)
axs[1].tick_params(axis='y', labelsize=16)
axs[1].text(klatWT[0]+90,1950,'WESTERN BOUNDARY',fontsize=12, weight='bold')
axs[1].text(klatET[3]+20,1950,'EASTERN BOUNDARY',fontsize=12, weight='bold')
axs[1].text(klonNT[0]+130,1950,'NORTHERN BOUNDARY',fontsize=12, weight='bold')
#colorbar:
cbar1=fig.colorbar(cax1,ax=axs[1],fraction=0.029, pad=0.02, ticks=np.arange(-0.6,0.7,0.1))
cbar1.ax.set_title('g/kg',size=16)
cbar1.outline.set_linewidth(0.4)
cbar1.ax.tick_params(labelsize=16,which='both')

cax2=axs[2].contourf(xxU,nc_Ubcl_p.depthu,U_anom_reshaped*1.e2,np.arange(-2.5,2.525,0.25),cmap='RdBu_r')
#axs[2].contour(xxV,nc_Vbcl_p.depthv,V_anom_reshaped*1.e-2,colors='k')
axs[2].pcolormesh(xxT,nc_Ubcl_p.depthu,msk_reshaped,vmin=-2,vmax=2,cmap='Greys',rasterized=True)
axs[2].plot([nWU, nWU],[depmax,0.],color='k',linewidth=1.9)
axs[2].plot([nWU+nNU, nWU+nNU],[depmax,0.],color='k',linewidth=1.9)
for kk in np.arange(np.size(klatWU)):
  axs[2].plot([klatWU[kk],klatWU[kk]],[depmax,0.],color='k',linestyle='--',linewidth=0.6)
  axs[2].plot([klatEU[kk],klatEU[kk]],[depmax,0.],color='k',linestyle='--',linewidth=0.6)
for kk in np.arange(np.size(klonNU)):
  axs[2].plot([klonNU[kk],klonNU[kk]],[depmax,0.],color='k',linestyle='--',linewidth=0.6)
axs[2].set_title('(c) Zonal velocity anomaly',fontsize=20,weight='bold')
axs[2].set_ylabel('Depth (m)',fontsize=18)
axs[2].invert_yaxis()
axs[2].set_ylim(depmax,0)
axs[2].set_xticks(np.append(np.append(klatWU,klonNU),np.flip(klatEU)))
axs[2].set_xticklabels(['75°S','70°S','65°S','60°S','135°W','120°W','105°W','90°W','60°S','65°S','70°S','75°S'])
axs[2].tick_params(axis='x', labelsize=16)
axs[2].tick_params(axis='y', labelsize=16)
axs[2].text(klatWU[0]+90,1950,'WESTERN BOUNDARY',fontsize=12, weight='bold')
axs[2].text(klatEU[3]+20,1950,'EASTERN BOUNDARY',fontsize=12, weight='bold')
axs[2].text(klonNU[0]+130,1950,'NORTHERN BOUNDARY',fontsize=12, weight='bold')
#colorbar:
cbar2=fig.colorbar(cax2,ax=axs[2],fraction=0.029, pad=0.02, ticks=np.arange(-2.5,2.55,0.5))
cbar2.ax.set_title('cm/s',size=16)
cbar2.outline.set_linewidth(0.4)
cbar2.ax.tick_params(labelsize=16,which='both')

#----------------------------------------------------------------------
fig.savefig('BDY_anom.jpg')
fig.savefig('BDY_anom.pdf')
