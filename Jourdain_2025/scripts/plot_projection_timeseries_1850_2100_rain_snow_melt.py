import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os

fig, axs = plt.subplots(nrows=2,ncols=1,figsize=(18.0,27.0))
axs = axs.ravel()

def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float,axis=0)
    ret[n:,:] = ret[n:,:] - ret[:-n,:]
    return ret[n - 1:,:] / n

def moving_average_1d(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

grd=xr.open_dataset('RCM_ice_regrid_04000m.nc2',decode_cf=False)
grd = grd.isel(x=slice(10,1510),y=slice(10,1510))
msk_isf = ( grd.ICE_MAR - grd.GROUND ) * grd.af2
msk_grd = grd.GROUND * grd.af2
print(msk_grd)

print('Total grounded ice sheet area = ',msk_grd.sum(dim=["x","y"]).values * 4. * 4. * 1.e-6,' million km2')
print('Total ice shelf area = ',msk_isf.sum(dim=["x","y"]).values * 4. * 4. * 1.e-6,' million km2')

# conversion to Gt/yr
#fac=4.e3*4.e3*1.e-12*86400*365.25
fac=4.e3*4.e3*1.e-12

ds=xr.open_dataset('../RCMs_summer_paper/in/MAR3.11-IPSL_ssp585-onbedmachinev2-4km_1980-2299.nc',drop_variables={'x','y'},decode_cf=False)

print(ds.SForg)
print(ds.SForg * msk_grd)
snow_grd = (ds.SForg * msk_grd * fac).sum(dim=["x","y"]).values
rain_grd = (ds.RForg * msk_grd * fac).sum(dim=["x","y"]).values
melt_grd = (ds.MEorg * msk_grd * fac).sum(dim=["x","y"]).values

snow_isf = (ds.SForg * msk_isf * fac).sum(dim=["x","y"]).values
rain_isf = (ds.RForg * msk_isf * fac).sum(dim=["x","y"]).values
melt_isf = (ds.MEorg * msk_isf * fac).sum(dim=["x","y"]).values

# running averages:
time_sm = moving_average_1d(np.arange(1980,2300),21)
snow_grd_sm = moving_average_1d(snow_grd,21)
rain_grd_sm = moving_average_1d(rain_grd,21)
melt_grd_sm = moving_average_1d(melt_grd,21)
snow_isf_sm = moving_average_1d(snow_isf,21)
rain_isf_sm = moving_average_1d(rain_isf,21)
melt_isf_sm = moving_average_1d(melt_isf,21)

# some diags:
print('##### grounded ice #####')
print(time_sm[10],' rain/precip =',rain_grd_sm[10]/(snow_grd_sm[10]+rain_grd_sm[10]),' rain/melt =',rain_grd_sm[10]/melt_grd_sm[10])
print(time_sm[110],' rain/precip =',rain_grd_sm[110]/(snow_grd_sm[110]+rain_grd_sm[110]),' rain/melt =',rain_grd_sm[110]/melt_grd_sm[110])
print(time_sm[210],' rain/precip =',rain_grd_sm[210]/(snow_grd_sm[210]+rain_grd_sm[210]),' rain/melt =',rain_grd_sm[210]/melt_grd_sm[210])
print('##### ice shelves #####')
print(time_sm[10],' rain/precip =',rain_isf_sm[10]/(snow_isf_sm[10]+rain_isf_sm[10]),' rain/melt =',rain_isf_sm[10]/melt_isf_sm[10])
print(time_sm[110],' rain/precip =',rain_isf_sm[110]/(snow_isf_sm[110]+rain_isf_sm[110]),' rain/melt =',rain_isf_sm[110]/melt_isf_sm[110])
print(time_sm[210],' rain/precip =',rain_isf_sm[210]/(snow_isf_sm[210]+rain_isf_sm[210]),' rain/melt =',rain_isf_sm[210]/melt_isf_sm[210])
print(time_sm[10],' snow =',snow_isf_sm[10],' rain =',rain_isf_sm[10],' melt=',melt_isf_sm[10])
print(time_sm[110],' snow =',snow_isf_sm[110],' rain =',rain_isf_sm[110],' melt=',melt_isf_sm[110])
print(time_sm[210],' snow =',snow_isf_sm[210],' rain =',rain_isf_sm[210],' melt=',melt_isf_sm[210])

# plots:
axs[0].plot(time_sm,snow_grd_sm,color='darkblue',linewidth=2.5,label='Snowfall')
axs[0].plot(time_sm,rain_grd_sm,color='orchid',linewidth=2.5,label='Rainfall')
axs[0].plot(time_sm,melt_grd_sm,color='orange',linewidth=2.5,label='Surface melting')

axs[1].plot(time_sm,snow_isf_sm,color='darkblue',linewidth=2.5,label='Snowfall')
axs[1].plot(time_sm,rain_isf_sm,color='orchid',linewidth=2.5,label='Rainfall')
axs[1].plot(time_sm,melt_isf_sm,color='orange',linewidth=2.5,label='Surface melting')

axs[0].set_title('(a) Grounded ice sheet',fontsize=28,fontweight='bold')
axs[1].set_title('(b) Ice shelves',fontsize=28,fontweight='bold')
axs[0].legend(loc='upper left',fontsize=26)
axs[0].set_ylim([0,4300])
axs[1].set_ylim([0,2990])
for kk in np.arange(2):
  axs[kk].set_ylabel('Gt/yr',fontsize=24)
  axs[kk].set_xlim([1990,2200])
  axs[kk].tick_params(axis='x', labelsize=24)
  axs[kk].tick_params(axis='y', labelsize=24)
  axs[kk].grid(color = 'black', linestyle = 'dotted', linewidth=1)

##########################################################################

fig.savefig("rain_snow_melt_IPSL_1980_2200.pdf")
