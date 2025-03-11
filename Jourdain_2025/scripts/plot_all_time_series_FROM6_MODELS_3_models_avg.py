import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

fig, axs = plt.subplots(nrows=3,ncols=2,figsize=(18.0,22.0))
axs = axs.ravel()

# input files
grd=xr.open_dataset('RCM_ice_regrid_04000m.nc2',decode_cf=False)
#-
model=['MPI-ESM1-2-HR',  'UKESM1-0-LL',  'CESM2']
member=[ 'r1i1p1f1', 'r1i1p1f2', 'r11i1p1f1' ]

for kmod in np.arange(3):
   
   exec("smhis_"+kmod.astype('str')+"=xr.open_dataset('MAR-'+model[kmod]+'_asmb_1980-2014_histo_regrid_04000m.nc',decode_cf=False)")
   exec("sm585_"+kmod.astype('str')+"=xr.open_dataset('MAR-'+model[kmod]+'_asmb_2015-2100_ssp585_regrid_04000m.nc',decode_cf=False)")
   exec("sm585EXT_"+kmod.astype('str')+"=xr.open_dataset('MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_2015-2100_ssp585_regrid_04000m_FROM_5_MODELS.nc',decode_cf=False)")
   exec("sm245_"+kmod.astype('str')+"=xr.open_dataset('MAR-'+model[kmod]+'_asmb_2015-2100_ssp245_regrid_04000m.nc',decode_cf=False)")
   exec("sm245EXT_"+kmod.astype('str')+"=xr.open_dataset('MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_2015-2100_ssp245_regrid_04000m_FROM_5_MODELS.nc',decode_cf=False)")
   exec("sm126_"+kmod.astype('str')+"=xr.open_dataset('MAR-'+model[kmod]+'_asmb_2015-2100_ssp126_regrid_04000m.nc',decode_cf=False)")
   exec("sm126EXT_"+kmod.astype('str')+"=xr.open_dataset('MAR-'+model[kmod]+'-'+member[kmod]+'_asmb_2015-2100_ssp126_regrid_04000m_FROM_5_MODELS.nc',decode_cf=False)")
   #-
   exec("mehis_"+kmod.astype('str')+"=xr.open_dataset('MAR-'+model[kmod]+'_ame_1980-2014_histo_regrid_04000m.nc',decode_cf=False)")
   exec("me585_"+kmod.astype('str')+"=xr.open_dataset('MAR-'+model[kmod]+'_ame_2015-2100_ssp585_regrid_04000m.nc',decode_cf=False)")
   exec("me585EXT_"+kmod.astype('str')+"=xr.open_dataset('MAR-'+model[kmod]+'-'+member[kmod]+'_ame_2015-2100_ssp585_regrid_04000m_FROM_5_MODELS.nc',decode_cf=False)")
   exec("me245_"+kmod.astype('str')+"=xr.open_dataset('MAR-'+model[kmod]+'_ame_2015-2100_ssp245_regrid_04000m.nc',decode_cf=False)")
   exec("me245EXT_"+kmod.astype('str')+"=xr.open_dataset('MAR-'+model[kmod]+'-'+member[kmod]+'_ame_2015-2100_ssp245_regrid_04000m_FROM_5_MODELS.nc',decode_cf=False)")
   exec("me126_"+kmod.astype('str')+"=xr.open_dataset('MAR-'+model[kmod]+'_ame_2015-2100_ssp126_regrid_04000m.nc',decode_cf=False)")
   exec("me126EXT_"+kmod.astype('str')+"=xr.open_dataset('MAR-'+model[kmod]+'-'+member[kmod]+'_ame_2015-2100_ssp126_regrid_04000m_FROM_5_MODELS.nc',decode_cf=False)")
   #-
   exec("ruhis_"+kmod.astype('str')+"=xr.open_dataset('MAR-'+model[kmod]+'_aru_1980-2014_histo_regrid_04000m.nc',decode_cf=False)")
   exec("ru585_"+kmod.astype('str')+"=xr.open_dataset('MAR-'+model[kmod]+'_aru_2015-2100_ssp585_regrid_04000m.nc',decode_cf=False)")
   exec("ru585EXT_"+kmod.astype('str')+"=xr.open_dataset('MAR-'+model[kmod]+'-'+member[kmod]+'_aru_2015-2100_ssp585_regrid_04000m_FROM_5_MODELS.nc',decode_cf=False)")
   exec("ru245_"+kmod.astype('str')+"=xr.open_dataset('MAR-'+model[kmod]+'_aru_2015-2100_ssp245_regrid_04000m.nc',decode_cf=False)")
   exec("ru245EXT_"+kmod.astype('str')+"=xr.open_dataset('MAR-'+model[kmod]+'-'+member[kmod]+'_aru_2015-2100_ssp245_regrid_04000m_FROM_5_MODELS.nc',decode_cf=False)")
   exec("ru126_"+kmod.astype('str')+"=xr.open_dataset('MAR-'+model[kmod]+'_aru_2015-2100_ssp126_regrid_04000m.nc',decode_cf=False)")
   exec("ru126EXT_"+kmod.astype('str')+"=xr.open_dataset('MAR-'+model[kmod]+'-'+member[kmod]+'_aru_2015-2100_ssp126_regrid_04000m_FROM_5_MODELS.nc',decode_cf=False)")

msk_isf = ( grd.ICE_MAR - grd.GROUND ) * grd.af2
msk_grn = grd.GROUND * grd.af2

# conversion to Gt/yr
fac=4.e3*4.e3*1.e-12*86400*365.25

#----- a -----
axs[0].plot([1980,2100],[0,0],linewidth=0.5,color='grey')
tmp1 = (smhis_0.asmb + smhis_1.asmb + smhis_2.asmb) * msk_grn * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[0].plot(np.arange(1980,2015),tmp2,color='black')
tmp1 = (sm585_0.asmb + sm585_1.asmb + sm585_2.asmb) * msk_grn * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[0].plot(np.arange(2015,2101),tmp2,color='black')
mean_2015_2100=np.mean(tmp2)
ymin, ymax = axs[0].get_ylim()
dy=ymax-ymin
axs[0].text(2020,ymax-0.05*dy,'SSP5-8.5 from RCMs  (mean='+np.array2string(mean_2015_2100,precision=1)+')',fontsize=14,color='black')
tmp1 = (sm585EXT_0.asmb + sm585EXT_1.asmb + sm585EXT_2.asmb) * msk_grn * fac / 3.; tmp3=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
bias_2015_2100=np.mean(tmp3-tmp2)
axs[0].text(2020,ymax-0.10*dy,'Emulated SSP5-8.5 (bias='+np.array2string(bias_2015_2100,precision=1)+')',fontsize=14,color='firebrick')
tmp1 = (sm585EXT_0.mmm_err + sm585EXT_1.mmm_err + sm585EXT_2.mmm_err) * msk_grn * fac * 1.96 / np.sqrt(5) / 3.; tmp4=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
#axs[0].fill_between(np.arange(2015,2101),tmp3-tmp4,tmp3+tmp4,color='firebrick',alpha=0.1)
axs[0].plot(np.arange(2015,2101),tmp3,color='firebrick')
#-
tmp1 = (sm245_0.asmb + sm245_1.asmb + sm245_2.asmb) * msk_grn * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[0].plot(np.arange(2015,2101),tmp2,color=[0.5,0.5,0.5])
mean_2015_2100=np.mean(tmp2)
axs[0].text(2020,ymax-0.15*dy,'SSP2-4.5 from RCMs  (mean='+np.array2string(mean_2015_2100,precision=1)+')',fontsize=14,color=[0.5,0.5,0.5])
tmp1 = (sm245EXT_0.asmb + sm245EXT_1.asmb + sm245EXT_2.asmb) * msk_grn * fac / 3.; tmp3=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
bias_2015_2100=np.mean(tmp3-tmp2)
axs[0].text(2020,ymax-0.20*dy,'Emulated SSP2-4.5 (bias='+np.array2string(bias_2015_2100,precision=1)+')',fontsize=14,color='red')
tmp1 = (sm245EXT_0.mmm_err + sm245EXT_1.mmm_err + sm245EXT_2.mmm_err) * msk_grn * fac * 1.96 / np.sqrt(5) / 3.; tmp4=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
#axs[0].fill_between(np.arange(2015,2101),tmp3-tmp4,tmp3+tmp4,color='red',alpha=0.1)
axs[0].plot(np.arange(2015,2101),tmp3,color='red')
#-
tmp1 = (sm126_0.asmb + sm126_1.asmb + sm126_2.asmb) * msk_grn * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[0].plot(np.arange(2015,2101),tmp2,color=[0.8,0.8,0.8])
mean_2015_2100=np.mean(tmp2)
axs[0].text(2020,ymax-0.25*dy,'SSP1-2.6 from RCMs  (mean='+np.array2string(mean_2015_2100,precision=1)+')',fontsize=14,color=[0.8,0.8,0.8])
tmp1 = (sm126EXT_0.asmb + sm126EXT_1.asmb + sm126EXT_2.asmb) * msk_grn * fac / 3.; tmp3=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
bias_2015_2100=np.mean(tmp3-tmp2)
axs[0].text(2020,ymax-0.30*dy,'Emulated SSP1-2.6 (bias='+np.array2string(bias_2015_2100,precision=1)+')',fontsize=14,color='orange')
tmp1 = (sm126EXT_0.mmm_err + sm126EXT_1.mmm_err + sm126EXT_2.mmm_err) * msk_grn * fac * 1.96 / np.sqrt(5) / 3.; tmp4=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
#axs[0].fill_between(np.arange(2015,2101),tmp3-tmp4,tmp3+tmp4,color='orange',alpha=0.1)
axs[0].plot(np.arange(2015,2101),tmp3,color='orange')
#-
axs[0].set_title('(a) SMB anomaly over the grounded ice sheet',fontsize=16,fontweight='bold')
axs[0].tick_params(axis='both', labelsize=12)
axs[0].set_ylabel('Gt/yr',fontsize=14)
axs[0].set_xlim([2015,2100])

#----- b -----
axs[1].plot([1980,2100],[0,0],linewidth=0.5,color='grey')
tmp1 = (smhis_0.asmb + smhis_1.asmb + smhis_2.asmb) * msk_isf * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[1].plot(np.arange(1980,2015),tmp2,color='black')
tmp1 = (sm585_0.asmb + sm585_1.asmb + sm585_2.asmb) * msk_isf * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[1].plot(np.arange(2015,2101),tmp2,color='black')
mean_2015_2100=np.mean(tmp2)
ymin, ymax = axs[1].get_ylim()
dy=ymax-ymin
axs[1].text(2020,ymin+0.60*dy,'SSP5-8.5 from RCMs  (mean='+np.array2string(mean_2015_2100,precision=1)+')',fontsize=14,color='black')
tmp1 = (sm585EXT_0.asmb + sm585EXT_1.asmb + sm585EXT_2.asmb) * msk_isf * fac / 3.; tmp3=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
bias_2015_2100=np.mean(tmp3-tmp2)
axs[1].text(2020,ymin+0.50*dy,'Emulated SSP5-8.5 (bias='+np.array2string(bias_2015_2100,precision=1)+')',fontsize=14,color='firebrick')
tmp1 = (sm585EXT_0.mmm_err + sm585EXT_1.mmm_err + sm585EXT_2.mmm_err) * msk_isf * fac * 1.96 / np.sqrt(5) / 3.; tmp4=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
#axs[1].fill_between(np.arange(2015,2101),tmp3-tmp4,tmp3+tmp4,color='firebrick',alpha=0.1)
axs[1].plot(np.arange(2015,2101),tmp3,color='firebrick')
#-
tmp1 = (sm245_0.asmb + sm245_1.asmb + sm245_2.asmb) * msk_isf * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[1].plot(np.arange(2015,2101),tmp2,color=[0.5,0.5,0.5])
mean_2015_2100=np.mean(tmp2)
axs[1].text(2020,ymin+0.40*dy,'SSP2-4.5 from RCMs  (mean='+np.array2string(mean_2015_2100,precision=1)+')',fontsize=14,color=[0.5,0.5,0.5])
tmp1 = (sm245EXT_0.asmb + sm245EXT_1.asmb + sm245EXT_2.asmb) * msk_isf * fac / 3.; tmp3=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
bias_2015_2100=np.mean(tmp3-tmp2)
axs[1].text(2020,ymin+0.30*dy,'Emulated SSP2-4.5 (bias='+np.array2string(bias_2015_2100,precision=1)+')',fontsize=14,color='red')
tmp1 = (sm245EXT_0.mmm_err + sm245EXT_1.mmm_err + sm245EXT_2.mmm_err) * msk_isf * fac * 1.96 / np.sqrt(5) / 3.; tmp4=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
#axs[1].fill_between(np.arange(2015,2101),tmp3-tmp4,tmp3+tmp4,color='red',alpha=0.1)
axs[1].plot(np.arange(2015,2101),tmp3,color='red')
#-
tmp1 = (sm126_0.asmb + sm126_1.asmb + sm126_2.asmb) * msk_isf * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[1].plot(np.arange(2015,2101),tmp2,color=[0.8,0.8,0.8])
mean_2015_2100=np.mean(tmp2)
axs[1].text(2020,ymin+0.20*dy,'SSP1-2.6 from RCMs  (mean='+np.array2string(mean_2015_2100,precision=1)+')',fontsize=14,color=[0.8,0.8,0.8])
tmp1 = (sm126EXT_0.asmb + sm126EXT_1.asmb + sm126EXT_2.asmb) * msk_isf * fac / 3.; tmp3=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
bias_2015_2100=np.mean(tmp3-tmp2)
axs[1].text(2020,ymin+0.10*dy,'Emulated SSP1-2.6 (bias='+np.array2string(bias_2015_2100,precision=1)+')',fontsize=14,color='orange')
tmp1 = (sm126EXT_0.mmm_err + sm126EXT_1.mmm_err + sm126EXT_2.mmm_err) * msk_isf * fac * 1.96 / np.sqrt(5) / 3.; tmp4=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
#axs[1].fill_between(np.arange(2015,2101),tmp3-tmp4,tmp3+tmp4,color='orange',alpha=0.1)
axs[1].plot(np.arange(2015,2101),tmp3,color='orange')
#-
axs[1].set_title('(b) SMB anomaly over the ice shelves',fontsize=16,fontweight='bold')
axs[1].tick_params(axis='both', labelsize=12)
axs[1].set_ylabel('Gt/yr',fontsize=14)
axs[1].set_xlim([2015,2100])

#----- c -----
axs[2].plot([1980,2100],[0,0],linewidth=0.5,color='grey')
tmp1 = (mehis_0.ame + mehis_1.ame + mehis_2.ame) * msk_grn * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[2].plot(np.arange(1980,2015),tmp2,color='black')
tmp1 = (me585_0.ame + me585_1.ame + me585_2.ame) * msk_grn * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[2].plot(np.arange(2015,2101),tmp2,color='black')
mean_2015_2100=np.mean(tmp2)
ymin, ymax = axs[2].get_ylim()
dy=ymax-ymin
axs[2].text(2020,ymax-0.07*dy,'SSP5-8.5 from RCMs  (mean='+np.array2string(mean_2015_2100,precision=1)+')',fontsize=14,color='black')
tmp1 = (me585EXT_0.ame + me585EXT_1.ame + me585EXT_2.ame) * msk_grn * fac / 3.; tmp3=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
bias_2015_2100=np.mean(tmp3-tmp2)
axs[2].text(2020,ymax-0.14*dy,'Emulated SSP5-8.5 (bias='+np.array2string(bias_2015_2100,precision=1)+')',fontsize=14,color='firebrick')
tmp1 = (me585EXT_0.mmm_err + me585EXT_1.mmm_err + me585EXT_2.mmm_err) * msk_grn * fac * 1.96 / np.sqrt(5) / 3.; tmp4=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
#axs[2].fill_between(np.arange(2015,2101),tmp3-tmp4,tmp3+tmp4,color='firebrick',alpha=0.1)
axs[2].plot(np.arange(2015,2101),tmp3,color='firebrick')
#-
tmp1 = (me245_0.ame + me245_1.ame + me245_2.ame) * msk_grn * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[2].plot(np.arange(2015,2101),tmp2,color=[0.5,0.5,0.5])
mean_2015_2100=np.mean(tmp2)
axs[2].text(2020,ymax-0.21*dy,'SSP2-4.5 from RCMs  (mean='+np.array2string(mean_2015_2100,precision=1)+')',fontsize=14,color=[0.5,0.5,0.5])
tmp1 = (me245EXT_0.ame + me245EXT_1.ame + me245EXT_2.ame) * msk_grn * fac / 3.; tmp3=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
bias_2015_2100=np.mean(tmp3-tmp2)
axs[2].text(2020,ymax-0.28*dy,'Emulated SSP2-4.5 (bias='+np.array2string(bias_2015_2100,precision=1)+')',fontsize=14,color='red')
tmp1 = (me245EXT_0.mmm_err + me245EXT_1.mmm_err + me245EXT_2.mmm_err) * msk_grn * fac * 1.96 / np.sqrt(5) / 3.; tmp4=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
#axs[2].fill_between(np.arange(2015,2101),tmp3-tmp4,tmp3+tmp4,color='red',alpha=0.1)
axs[2].plot(np.arange(2015,2101),tmp3,color='red')
#-
tmp1 = (me126_0.ame + me126_1.ame + me126_2.ame) * msk_grn * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[2].plot(np.arange(2015,2101),tmp2,color=[0.8,0.8,0.8])
mean_2015_2100=np.mean(tmp2)
axs[2].text(2020,ymax-0.35*dy,'SSP1-2.6 from RCMs  (mean='+np.array2string(mean_2015_2100,precision=1)+')',fontsize=14,color=[0.8,0.8,0.8])
tmp1 = (me126EXT_0.ame + me126EXT_1.ame + me126EXT_2.ame) * msk_grn * fac / 3.; tmp3=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
bias_2015_2100=np.mean(tmp3-tmp2)
axs[2].text(2020,ymax-0.42*dy,'Emulated SSP1-2.6 (bias='+np.array2string(bias_2015_2100,precision=1)+')',fontsize=14,color='orange')
tmp1 = (me126EXT_0.mmm_err + me126EXT_1.mmm_err + me126EXT_2.mmm_err) * msk_grn * fac * 1.96 / np.sqrt(5) / 3.; tmp4=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
#axs[2].fill_between(np.arange(2015,2101),tmp3-tmp4,tmp3+tmp4,color='orange',alpha=0.1)
axs[2].plot(np.arange(2015,2101),tmp3,color='orange')
#-
axs[2].set_title('(c) Melting anomaly over the grounded ice sheet',fontsize=16,fontweight='bold')
axs[2].tick_params(axis='both', labelsize=12)
axs[2].set_ylabel('Gt/yr',fontsize=14)
axs[2].set_xlim([2015,2100])

#----- d -----
axs[3].plot([1980,2100],[0,0],linewidth=0.5,color='grey')
tmp1 = (mehis_0.ame + mehis_1.ame + mehis_2.ame) * msk_isf * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[3].plot(np.arange(1980,2015),tmp2,color='black')
tmp1 = (me585_0.ame + me585_1.ame + me585_2.ame) * msk_isf * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[3].plot(np.arange(2015,2101),tmp2,color='black')
mean_2015_2100=np.mean(tmp2)
ymin, ymax = axs[3].get_ylim()
dy=ymax-ymin
axs[3].text(2020,ymax-0.07*dy,'SSP5-8.5 from RCMs  (mean='+np.array2string(mean_2015_2100,precision=1)+')',fontsize=14,color='black')
tmp1 = (me585EXT_0.ame + me585EXT_1.ame + me585EXT_2.ame) * msk_isf * fac / 3.; tmp3=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
bias_2015_2100=np.mean(tmp3-tmp2)
axs[3].text(2020,ymax-0.14*dy,'Emulated SSP5-8.5 (bias='+np.array2string(bias_2015_2100,precision=1)+')',fontsize=14,color='firebrick')
tmp1 = (me585EXT_0.mmm_err + me585EXT_1.mmm_err + me585EXT_2.mmm_err) * msk_isf * fac * 1.96 / np.sqrt(5) / 3.; tmp4=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
#axs[3].fill_between(np.arange(2015,2101),tmp3-tmp4,tmp3+tmp4,color='firebrick',alpha=0.1)
axs[3].plot(np.arange(2015,2101),tmp3,color='firebrick')
#-
tmp1 = (me245_0.ame + me245_1.ame + me245_2.ame) * msk_isf * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[3].plot(np.arange(2015,2101),tmp2,color=[0.5,0.5,0.5])
mean_2015_2100=np.mean(tmp2)
axs[3].text(2020,ymax-0.21*dy,'SSP2-4.5 from RCMs  (mean='+np.array2string(mean_2015_2100,precision=1)+')',fontsize=14,color=[0.5,0.5,0.5])
tmp1 = (me245EXT_0.ame + me245EXT_1.ame + me245EXT_2.ame) * msk_isf * fac / 3.; tmp3=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
bias_2015_2100=np.mean(tmp3-tmp2)
axs[3].text(2020,ymax-0.28*dy,'Emulated SSP2-4.5 (bias='+np.array2string(bias_2015_2100,precision=1)+')',fontsize=14,color='red')
tmp1 = (me245EXT_0.mmm_err + me245EXT_1.mmm_err + me245EXT_2.mmm_err) * msk_isf * fac * 1.96 / np.sqrt(5) / 3.; tmp4=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
#axs[3].fill_between(np.arange(2015,2101),tmp3-tmp4,tmp3+tmp4,color='red',alpha=0.1)
axs[3].plot(np.arange(2015,2101),tmp3,color='red')
#-
tmp1 = (me126_0.ame + me126_1.ame + me126_2.ame) * msk_isf * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[3].plot(np.arange(2015,2101),tmp2,color=[0.8,0.8,0.8])
mean_2015_2100=np.mean(tmp2)
axs[3].text(2020,ymax-0.35*dy,'SSP1-2.6 from RCMs  (mean='+np.array2string(mean_2015_2100,precision=1)+')',fontsize=14,color=[0.8,0.8,0.8])
tmp1 = (me126EXT_0.ame + me126EXT_1.ame + me126EXT_2.ame) * msk_isf * fac / 3.; tmp3=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
bias_2015_2100=np.mean(tmp3-tmp2)
axs[3].text(2020,ymax-0.42*dy,'Emulated SSP1-2.6 (bias='+np.array2string(bias_2015_2100,precision=1)+')',fontsize=14,color='orange')
tmp1 = (me126EXT_0.mmm_err + me126EXT_1.mmm_err + me126EXT_2.mmm_err) * msk_isf * fac * 1.96 / np.sqrt(5) / 3.; tmp4=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
#axs[3].fill_between(np.arange(2015,2101),tmp3-tmp4,tmp3+tmp4,color='orange',alpha=0.1)
axs[3].plot(np.arange(2015,2101),tmp3,color='orange')
#-
axs[3].set_title('(d) Melting anomaly over the ice shelves',fontsize=16,fontweight='bold')
axs[3].tick_params(axis='both', labelsize=12)
axs[3].set_ylabel('Gt/yr',fontsize=14)
axs[3].set_xlim([2015,2100])

#----- e -----
axs[4].plot([1980,2100],[0,0],linewidth=0.5,color='grey')
tmp1 = (ruhis_0.aru + ruhis_1.aru + ruhis_2.aru) * msk_grn * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[4].plot(np.arange(1980,2015),tmp2,color='black')
tmp1 = (ru585_0.aru + ru585_1.aru + ru585_2.aru) * msk_grn * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[4].plot(np.arange(2015,2101),tmp2,color='black')
mean_2015_2100=np.mean(tmp2)
ymin, ymax = axs[4].get_ylim()
dy=ymax-ymin
axs[4].text(2020,ymax-0.07*dy,'SSP5-8.5 from RCMs  (mean='+np.array2string(mean_2015_2100,precision=1)+')',fontsize=14,color='black')
tmp1 = (ru585EXT_0.aru + ru585EXT_1.aru + ru585EXT_2.aru) * msk_grn * fac / 3.; tmp3=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
bias_2015_2100=np.mean(tmp3-tmp2)
axs[4].text(2020,ymax-0.14*dy,'Emulated SSP5-8.5 (bias='+np.array2string(bias_2015_2100,precision=1)+')',fontsize=14,color='firebrick')
tmp1 = (ru585EXT_0.mmm_err + ru585EXT_1.mmm_err + ru585EXT_2.mmm_err) * msk_grn * fac * 1.96 / np.sqrt(5) / 3.; tmp4=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
#axs[4].fill_between(np.arange(2015,2101),tmp3-tmp4,tmp3+tmp4,color='firebrick',alpha=0.1)
axs[4].plot(np.arange(2015,2101),tmp3,color='firebrick')
#-
tmp1 = (ru245_0.aru + ru245_1.aru + ru245_2.aru) * msk_grn * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[4].plot(np.arange(2015,2101),tmp2,color=[0.5,0.5,0.5])
mean_2015_2100=np.mean(tmp2)
axs[4].text(2020,ymax-0.21*dy,'SSP2-4.5 from RCMs  (mean='+np.array2string(mean_2015_2100,precision=1)+')',fontsize=14,color=[0.5,0.5,0.5])
tmp1 = (ru245EXT_0.aru + ru245EXT_1.aru + ru245EXT_2.aru) * msk_grn * fac / 3.; tmp3=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
bias_2015_2100=np.mean(tmp3-tmp2)
axs[4].text(2020,ymax-0.28*dy,'Emulated SSP2-4.5 (bias='+np.array2string(bias_2015_2100,precision=1)+')',fontsize=14,color='red')
tmp1 = (ru245EXT_0.mmm_err + ru245EXT_1.mmm_err + ru245EXT_2.mmm_err) * msk_grn * fac * 1.96 / np.sqrt(5) / 3.; tmp4=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
#axs[4].fill_between(np.arange(2015,2101),tmp3-tmp4,tmp3+tmp4,color='red',alpha=0.1)
axs[4].plot(np.arange(2015,2101),tmp3,color='red')
#-
tmp1 = (ru126_0.aru + ru126_1.aru + ru126_2.aru) * msk_grn * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[4].plot(np.arange(2015,2101),tmp2,color=[0.8,0.8,0.8])
mean_2015_2100=np.mean(tmp2)
axs[4].text(2020,ymax-0.35*dy,'SSP1-2.6 from RCMs  (mean='+np.array2string(mean_2015_2100,precision=1)+')',fontsize=14,color=[0.8,0.8,0.8])
tmp1 = (ru126EXT_0.aru + ru126EXT_1.aru + ru126EXT_2.aru) * msk_grn * fac / 3.; tmp3=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
bias_2015_2100=np.mean(tmp3-tmp2)
axs[4].text(2020,ymax-0.42*dy,'Emulated SSP1-2.6 (bias='+np.array2string(bias_2015_2100,precision=1)+')',fontsize=14,color='orange')
tmp1 = (ru126EXT_0.mmm_err + ru126EXT_1.mmm_err + ru126EXT_2.mmm_err) * msk_grn * fac * 1.96 / np.sqrt(5) / 3.; tmp4=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
#axs[4].fill_between(np.arange(2015,2101),tmp3-tmp4,tmp3+tmp4,color='orange',alpha=0.1)
axs[4].plot(np.arange(2015,2101),tmp3,color='orange')
#-
axs[4].set_title('(e) Anom. excess liquid water over the grounded i. s.',fontsize=16,fontweight='bold')
axs[4].tick_params(axis='both', labelsize=12)
axs[4].set_ylabel('Gt/yr',fontsize=14)
axs[4].set_xlim([2015,2100])

#----- f -----
axs[5].plot([1980,2100],[0,0],linewidth=0.5,color='grey')
tmp1 = (ruhis_0.aru + ruhis_1.aru + ruhis_2.aru) * msk_isf * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[5].plot(np.arange(1980,2015),tmp2,color='black')
tmp1 = (ru585_0.aru + ru585_1.aru + ru585_2.aru) * msk_isf * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[5].plot(np.arange(2015,2101),tmp2,color='black')
mean_2015_2100=np.mean(tmp2)
ymin, ymax = axs[5].get_ylim()
dy=ymax-ymin
axs[5].text(2020,ymax-0.07*dy,'SSP5-8.5 from RCMs  (mean='+np.array2string(mean_2015_2100,precision=1)+')',fontsize=14,color='black')
tmp1 = (ru585EXT_0.aru + ru585EXT_1.aru + ru585EXT_2.aru) * msk_isf * fac / 3.; tmp3=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
bias_2015_2100=np.mean(tmp3-tmp2)
axs[5].text(2020,ymax-0.14*dy,'Emulated SSP5-8.5 (bias='+np.array2string(bias_2015_2100,precision=1)+')',fontsize=14,color='firebrick')
tmp1 = (ru585EXT_0.mmm_err + ru585EXT_1.mmm_err + ru585EXT_2.mmm_err) * msk_isf * fac * 1.96 / np.sqrt(5) / 3.; tmp4=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
#axs[5].fill_between(np.arange(2015,2101),tmp3-tmp4,tmp3+tmp4,color='firebrick',alpha=0.1)
axs[5].plot(np.arange(2015,2101),tmp3,color='firebrick')
#-
tmp1 = (ru245_0.aru + ru245_1.aru + ru245_2.aru) * msk_isf * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[5].plot(np.arange(2015,2101),tmp2,color=[0.5,0.5,0.5])
mean_2015_2100=np.mean(tmp2)
axs[5].text(2020,ymax-0.21*dy,'SSP2-4.5 from RCMs  (mean='+np.array2string(mean_2015_2100,precision=1)+')',fontsize=14,color=[0.5,0.5,0.5])
tmp1 = (ru245EXT_0.aru + ru245EXT_1.aru + ru245EXT_2.aru) * msk_isf * fac / 3.; tmp3=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
bias_2015_2100=np.mean(tmp3-tmp2)
axs[5].text(2020,ymax-0.28*dy,'Emulated SSP2-4.5 (bias='+np.array2string(bias_2015_2100,precision=1)+')',fontsize=14,color='red')
tmp1 = (ru245EXT_0.mmm_err + ru245EXT_1.mmm_err + ru245EXT_2.mmm_err) * msk_isf * fac * 1.96 / np.sqrt(5) / 3.; tmp4=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
#axs[5].fill_between(np.arange(2015,2101),tmp3-tmp4,tmp3+tmp4,color='red',alpha=0.1)
axs[5].plot(np.arange(2015,2101),tmp3,color='red')
#-
tmp1 = (ru126_0.aru + ru126_1.aru + ru126_2.aru) * msk_isf * fac / 3.; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[5].plot(np.arange(2015,2101),tmp2,color=[0.8,0.8,0.8])
mean_2015_2100=np.mean(tmp2)
axs[5].text(2020,ymax-0.35*dy,'SSP1-2.6 from RCMs  (mean='+np.array2string(mean_2015_2100,precision=1)+')',fontsize=14,color=[0.8,0.8,0.8])
tmp1 = (ru126EXT_0.aru + ru126EXT_1.aru + ru126EXT_2.aru) * msk_isf * fac / 3.; tmp3=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
bias_2015_2100=np.mean(tmp3-tmp2)
axs[5].text(2020,ymax-0.42*dy,'Emulated SSP1-2.6 (bias='+np.array2string(bias_2015_2100,precision=1)+')',fontsize=14,color='orange')
tmp1 = (ru126EXT_0.mmm_err + ru126EXT_1.mmm_err + ru126EXT_2.mmm_err) * msk_isf * fac * 1.96 / np.sqrt(5) / 3.; tmp4=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
#axs[5].fill_between(np.arange(2015,2101),tmp3-tmp4,tmp3+tmp4,color='orange',alpha=0.1)
axs[5].plot(np.arange(2015,2101),tmp3,color='orange')
#-
axs[5].set_title('(f) Anom. excess liquid water over the ice shelves',fontsize=16,fontweight='bold')
axs[5].tick_params(axis='both', labelsize=12)
axs[5].set_ylabel('Gt/yr',fontsize=14)
axs[5].set_xlim([2015,2100])

#-----

figname='timeseries_CESM2_UKESM_MPI_FROM_5_MODELS_new.pdf'

fig.savefig(figname)
