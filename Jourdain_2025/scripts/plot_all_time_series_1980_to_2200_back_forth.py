import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

fig, axs = plt.subplots(nrows=3,ncols=2,figsize=(18.0,22.0))
axs = axs.ravel()

# colors (RGB)
col050=np.array([255,190,0])/255
col090=np.array([100,0,180])/255
col060=0.75*col050+0.25*col090
col070=0.50*col050+0.50*col090
col080=0.25*col050+0.75*col090

# input files
grd=xr.open_dataset('RCM_ice_regrid_04000m.nc2',decode_cf=False)
smhis=xr.open_dataset('MAR-IPSL-CM6A-LR_asmb_1980-2014_histo_regrid_04000m.nc',decode_cf=False)
sm585=xr.open_dataset('MAR-IPSL-CM6A-LR_asmb_2015-2200_ssp585_regrid_04000m.nc',decode_cf=False)
smEXT050=xr.open_dataset('MAR-IPSL-CM6A-LR_asmb_2101-2300_ssp585_regrid_04000m_EXTENDED_MOA050.nc',decode_cf=False)
smEXT060=xr.open_dataset('MAR-IPSL-CM6A-LR_asmb_2101-2300_ssp585_regrid_04000m_EXTENDED_MOA060.nc',decode_cf=False)
smEXT070=xr.open_dataset('MAR-IPSL-CM6A-LR_asmb_2101-2300_ssp585_regrid_04000m_EXTENDED_MOA070.nc',decode_cf=False)
smEXT080=xr.open_dataset('MAR-IPSL-CM6A-LR_asmb_2101-2300_ssp585_regrid_04000m_EXTENDED_MOA080.nc',decode_cf=False)
smEXT090=xr.open_dataset('MAR-IPSL-CM6A-LR_asmb_2101-2300_ssp585_regrid_04000m_EXTENDED_MOA090.nc',decode_cf=False)
smBAC050=xr.open_dataset('MAR-IPSL-CM6A-LR_asmb_2015-2100_ssp585_regrid_04000m_EXTENDED_BACK_MOA050.nc',decode_cf=False)
smBAC060=xr.open_dataset('MAR-IPSL-CM6A-LR_asmb_2015-2100_ssp585_regrid_04000m_EXTENDED_BACK_MOA060.nc',decode_cf=False)
smBAC070=xr.open_dataset('MAR-IPSL-CM6A-LR_asmb_2015-2100_ssp585_regrid_04000m_EXTENDED_BACK_MOA070.nc',decode_cf=False)
smBAC080=xr.open_dataset('MAR-IPSL-CM6A-LR_asmb_2015-2100_ssp585_regrid_04000m_EXTENDED_BACK_MOA080.nc',decode_cf=False)
smBAC090=xr.open_dataset('MAR-IPSL-CM6A-LR_asmb_2015-2100_ssp585_regrid_04000m_EXTENDED_BACK_MOA090.nc',decode_cf=False)
mehis=xr.open_dataset('MAR-IPSL-CM6A-LR_ame_1980-2014_histo_regrid_04000m.nc',decode_cf=False)
me585=xr.open_dataset('MAR-IPSL-CM6A-LR_ame_2015-2200_ssp585_regrid_04000m.nc',decode_cf=False)
meEXT050=xr.open_dataset('MAR-IPSL-CM6A-LR_ame_2101-2300_ssp585_regrid_04000m_EXTENDED_MOA050.nc',decode_cf=False)
meEXT060=xr.open_dataset('MAR-IPSL-CM6A-LR_ame_2101-2300_ssp585_regrid_04000m_EXTENDED_MOA060.nc',decode_cf=False)
meEXT070=xr.open_dataset('MAR-IPSL-CM6A-LR_ame_2101-2300_ssp585_regrid_04000m_EXTENDED_MOA070.nc',decode_cf=False)
meEXT080=xr.open_dataset('MAR-IPSL-CM6A-LR_ame_2101-2300_ssp585_regrid_04000m_EXTENDED_MOA080.nc',decode_cf=False)
meEXT090=xr.open_dataset('MAR-IPSL-CM6A-LR_ame_2101-2300_ssp585_regrid_04000m_EXTENDED_MOA090.nc',decode_cf=False)
meBAC050=xr.open_dataset('MAR-IPSL-CM6A-LR_ame_2015-2100_ssp585_regrid_04000m_EXTENDED_BACK_MOA050.nc',decode_cf=False)
meBAC060=xr.open_dataset('MAR-IPSL-CM6A-LR_ame_2015-2100_ssp585_regrid_04000m_EXTENDED_BACK_MOA060.nc',decode_cf=False)
meBAC070=xr.open_dataset('MAR-IPSL-CM6A-LR_ame_2015-2100_ssp585_regrid_04000m_EXTENDED_BACK_MOA070.nc',decode_cf=False)
meBAC080=xr.open_dataset('MAR-IPSL-CM6A-LR_ame_2015-2100_ssp585_regrid_04000m_EXTENDED_BACK_MOA080.nc',decode_cf=False)
meBAC090=xr.open_dataset('MAR-IPSL-CM6A-LR_ame_2015-2100_ssp585_regrid_04000m_EXTENDED_BACK_MOA090.nc',decode_cf=False)
ruhis=xr.open_dataset('MAR-IPSL-CM6A-LR_aru_1980-2014_histo_regrid_04000m.nc',decode_cf=False)
ru585=xr.open_dataset('MAR-IPSL-CM6A-LR_aru_2015-2200_ssp585_regrid_04000m.nc',decode_cf=False)
ruEXT050=xr.open_dataset('MAR-IPSL-CM6A-LR_aru_2101-2300_ssp585_regrid_04000m_EXTENDED_MOA050.nc',decode_cf=False)
ruEXT060=xr.open_dataset('MAR-IPSL-CM6A-LR_aru_2101-2300_ssp585_regrid_04000m_EXTENDED_MOA060.nc',decode_cf=False)
ruEXT070=xr.open_dataset('MAR-IPSL-CM6A-LR_aru_2101-2300_ssp585_regrid_04000m_EXTENDED_MOA070.nc',decode_cf=False)
ruEXT080=xr.open_dataset('MAR-IPSL-CM6A-LR_aru_2101-2300_ssp585_regrid_04000m_EXTENDED_MOA080.nc',decode_cf=False)
ruEXT090=xr.open_dataset('MAR-IPSL-CM6A-LR_aru_2101-2300_ssp585_regrid_04000m_EXTENDED_MOA090.nc',decode_cf=False)
ruBAC050=xr.open_dataset('MAR-IPSL-CM6A-LR_aru_2015-2100_ssp585_regrid_04000m_EXTENDED_BACK_MOA050.nc',decode_cf=False)
ruBAC060=xr.open_dataset('MAR-IPSL-CM6A-LR_aru_2015-2100_ssp585_regrid_04000m_EXTENDED_BACK_MOA060.nc',decode_cf=False)
ruBAC070=xr.open_dataset('MAR-IPSL-CM6A-LR_aru_2015-2100_ssp585_regrid_04000m_EXTENDED_BACK_MOA070.nc',decode_cf=False)
ruBAC080=xr.open_dataset('MAR-IPSL-CM6A-LR_aru_2015-2100_ssp585_regrid_04000m_EXTENDED_BACK_MOA080.nc',decode_cf=False)
ruBAC090=xr.open_dataset('MAR-IPSL-CM6A-LR_aru_2015-2100_ssp585_regrid_04000m_EXTENDED_BACK_MOA090.nc',decode_cf=False)

# IPSL CMIP outputs:
IPSL_CMIP_cli=xr.open_dataset('/data/njourdain/DATA_PROTECT/SMB/smb_Clim_IPSL-CM6A-LR_historical_r1i1p1f1_1995_2014.nc',decode_cf=False) 
IPSL_CMIP_ssp=xr.open_dataset('/data/njourdain/DATA_PROTECT/SMB/smb_Ayr_IPSL-CM6A-LR_ssp585_r1i1p1f1_2015_2300.nc',decode_cf=False) 

msk_isf = ( grd.ICE_MAR - grd.GROUND ) * grd.af2
msk_grn = grd.GROUND * grd.af2

# conversion to Gt/yr
fac=4.e3*4.e3*1.e-12*86400*365.25

#----- a -----
axs[0].fill([2080.5,2100.5,2100.5,2080.5,2080.5],[-500,-500,1550,1550,-500],color='lightgrey',edgecolor=None)
axs[0].plot([1980,2200],[0,0],linewidth=0.5,color='k')
xrp1 = (IPSL_CMIP_ssp.smb.isel(time=slice(0,185)) - IPSL_CMIP_cli.smb.squeeze(drop=True)) * msk_grn * fac
xrp2 = xrp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[0].plot(np.arange(2015,2200),xrp2,color='black',linestyle='--',linewidth=0.8)
tmp1 = smhis.asmb * msk_grn * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[0].plot(np.arange(1980,2015),tmp2,color='black')
tmp1 = sm585.asmb * msk_grn * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[0].plot(np.arange(2015,2200),tmp2,color='black')
#-
tmp1 = smEXT050.asmb * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[0].plot(np.arange(2101,2200),tmp3,color=col050)
bias_2101_2120_1=np.mean(tmp3[0:20]-tmp2[86:106])
#-
tmp1 = smEXT060.asmb * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[0].plot(np.arange(2101,2200),tmp3,color=col060)
bias_2101_2120_2=np.mean(tmp3[0:20]-tmp2[86:106])
#-
tmp1 = smEXT070.asmb * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[0].plot(np.arange(2101,2200),tmp3,color=col070)
bias_2101_2120_3=np.mean(tmp3[0:20]-tmp2[86:106])
#-
tmp1 = smEXT080.asmb * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[0].plot(np.arange(2101,2200),tmp3,color=col080)
bias_2101_2120_4=np.mean(tmp3[0:20]-tmp2[86:106])
#-
tmp1 = smEXT090.asmb * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[0].plot(np.arange(2101,2200),tmp3,color=col090)
bias_2101_2120_5=np.mean(tmp3[0:20]-tmp2[86:106])
#-----
tmp1 = smBAC050.asmb * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[0].plot(np.arange(2015,2081),tmp3,color=col050)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[0].text(2082,225,'r = 0.5  ('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120_1,precision=1)+')',fontsize=14,color=col050)
#-
tmp1 = smBAC060.asmb * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[0].plot(np.arange(2015,2081),tmp3,color=col060)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[0].text(2082,75,'r = 0.6  ('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120_2,precision=1)+')',fontsize=14,color=col060)
#-
tmp1 = smBAC070.asmb * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[0].plot(np.arange(2015,2081),tmp3,color=col070)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[0].text(2082,-75,'r = 0.7  ('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120_3,precision=1)+')',fontsize=14,color=col070)
#-
tmp1 = smBAC080.asmb * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[0].plot(np.arange(2015,2081),tmp3,color=col080)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[0].text(2082,-225,'r = 0.8  ('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120_4,precision=1)+')',fontsize=14,color=col080)
#-
tmp1 = smBAC090.asmb * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[0].plot(np.arange(2015,2081),tmp3,color=col090)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[0].text(2082,-375,'r = 0.9  ('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120_5,precision=1)+')',fontsize=14,color=col090)
#-----
axs[0].set_title('(a) SMB anomaly over the grounded ice sheet',fontsize=16,weight='bold')
axs[0].tick_params(axis='both', labelsize=12)
axs[0].set_ylabel('Gt/yr',fontsize=14)
axs[0].set_xlim([2015,2200])
axs[0].set_ylim([-500,1400])

#----- b -----
axs[1].fill([2080.5,2100.5,2100.5,2080.5,2080.5],[-4000,-4000,400,400,-4000],color='lightgrey',edgecolor=None)
axs[1].plot([1980,2200],[0,0],linewidth=0.5,color='k')
xrp1 = (IPSL_CMIP_ssp.smb.isel(time=slice(0,185)) - IPSL_CMIP_cli.smb.squeeze(drop=True)) * msk_isf * fac
xrp2 = xrp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[1].plot(np.arange(2015,2200),xrp2,color='black',linestyle='--',linewidth=0.8)
tmp1 = smhis.asmb * msk_isf * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[1].plot(np.arange(1980,2015),tmp2,color='black')
tmp1 = sm585.asmb * msk_isf * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[1].plot(np.arange(2015,2200),tmp2,color='black')
#-
tmp1 = smEXT050.asmb * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[1].plot(np.arange(2101,2200),tmp3,color=col050)
bias_2101_2120_1=np.mean(tmp3[0:20]-tmp2[86:106])
#-
tmp1 = smEXT060.asmb * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[1].plot(np.arange(2101,2200),tmp3,color=col060)
bias_2101_2120_2=np.mean(tmp3[0:20]-tmp2[86:106])
#-
tmp1 = smEXT070.asmb * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[1].plot(np.arange(2101,2200),tmp3,color=col070)
bias_2101_2120_3=np.mean(tmp3[0:20]-tmp2[86:106])
#-
tmp1 = smEXT080.asmb * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[1].plot(np.arange(2101,2200),tmp3,color=col080)
bias_2101_2120_4=np.mean(tmp3[0:20]-tmp2[86:106])
#-
tmp1 = smEXT090.asmb * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[1].plot(np.arange(2101,2200),tmp3,color=col090)
bias_2101_2120_5=np.mean(tmp3[0:20]-tmp2[86:106])
#-----
tmp1 = smBAC050.asmb * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[1].plot(np.arange(2015,2081),tmp3,color=col050)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[1].text(2020,-1000,'r = 0.5  ('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120_1,precision=1)+')',fontsize=14,color=col050)
#-
tmp1 = smBAC060.asmb * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[1].plot(np.arange(2015,2081),tmp3,color=col060)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[1].text(2020,-1500,'r = 0.6  ('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120_2,precision=1)+')',fontsize=14,color=col060)
#-
tmp1 = smBAC070.asmb * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[1].plot(np.arange(2015,2081),tmp3,color=col070)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[1].text(2020,-2000,'r = 0.7  ('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120_3,precision=1)+')',fontsize=14,color=col070)
#-
tmp1 = smBAC080.asmb * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[1].plot(np.arange(2015,2081),tmp3,color=col080)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[1].text(2020,-2500,'r = 0.8  ('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120_4,precision=1)+')',fontsize=14,color=col080)
#-
tmp1 = smBAC090.asmb * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[1].plot(np.arange(2015,2081),tmp3,color=col090)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[1].text(2020,-3000,'r = 0.9  ('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120_5,precision=1)+')',fontsize=14,color=col090)
#-----
axs[1].set_title('(b) SMB anomaly over the ice shelves',fontsize=16,weight='bold')
axs[1].tick_params(axis='both', labelsize=12)
axs[1].set_ylabel('Gt/yr',fontsize=14)
axs[1].set_xlim([2015,2200])
axs[1].set_ylim([-4000,400])

#----- c -----
axs[2].fill([2080.5,2100.5,2100.5,2080.5,2080.5],[-100,-100,3900,3900,-100],color='lightgrey',edgecolor=None)
axs[2].plot([1980,2200],[0,0],linewidth=0.5,color='k')
tmp1 = mehis.ame * msk_grn * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[2].plot(np.arange(1980,2015),tmp2,color='black')
tmp1 = me585.ame * msk_grn * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[2].plot(np.arange(2015,2200),tmp2,color='black')
#-
tmp1 = meEXT070.ame * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[2].plot(np.arange(2101,2200),tmp3,color=col070)
bias_2101_2120=np.mean(tmp3[0:20]-tmp2[86:106])
#-----
tmp1 = meBAC070.ame * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[2].plot(np.arange(2015,2081),tmp3,color=col070)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[2].text(2030,2500,'('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120,precision=1)+')',fontsize=14,color=col070)
#-----
axs[2].set_title('(c) Melting anomaly over the grounded ice sheet',fontsize=16,weight='bold')
axs[2].tick_params(axis='both', labelsize=12)
axs[2].set_ylabel('Gt/yr',fontsize=14)
axs[2].set_xlim([2015,2200])
axs[2].set_ylim([-100,3900])

#----- d -----
axs[3].fill([2080.5,2100.5,2100.5,2080.5,2080.5],[-200,-200,5100,5100,-200],color='lightgrey',edgecolor=None)
axs[3].plot([1980,2200],[0,0],linewidth=0.5,color='k')
tmp1 = mehis.ame * msk_isf * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[3].plot(np.arange(1980,2015),tmp2,color='black')
tmp1 = me585.ame * msk_isf * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[3].plot(np.arange(2015,2200),tmp2,color='black')
#-
tmp1 = meEXT070.ame * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[3].plot(np.arange(2101,2200),tmp3,color=col070)
bias_2101_2120=np.mean(tmp3[0:20]-tmp2[86:106])
#-----
tmp1 = meBAC070.ame * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[3].plot(np.arange(2015,2081),tmp3,color=col070)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[3].text(2030,3400,'('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120,precision=1)+')',fontsize=14,color=col070)
#-----
axs[3].set_title('(d) Melting anomaly over the ice shelves',fontsize=16,weight='bold')
axs[3].tick_params(axis='both', labelsize=12)
axs[3].set_ylabel('Gt/yr',fontsize=14)
axs[3].set_xlim([2015,2200])
axs[3].set_ylim([-200,5100])

#----- e -----
axs[4].fill([2080.5,2100.5,2100.5,2080.5,2080.5],[-100,-100,3400,3400,-100],color='lightgrey',edgecolor=None)
axs[4].plot([1980,2200],[0,0],linewidth=0.5,color='k')
tmp1 = ruhis.aru * msk_grn * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[4].plot(np.arange(1980,2015),tmp2,color='black')
tmp1 = ru585.aru * msk_grn * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[4].plot(np.arange(2015,2200),tmp2,color='black')
#-
tmp1 = ruEXT050.aru * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[4].plot(np.arange(2101,2200),tmp3,color=col050)
bias_2101_2120_1=np.mean(tmp3[0:20]-tmp2[86:106])
#-
tmp1 = ruEXT060.aru * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[4].plot(np.arange(2101,2200),tmp3,color=col060)
bias_2101_2120_2=np.mean(tmp3[0:20]-tmp2[86:106])
#-
tmp1 = ruEXT070.aru * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[4].plot(np.arange(2101,2200),tmp3,color=col070)
bias_2101_2120_3=np.mean(tmp3[0:20]-tmp2[86:106])
#-
tmp1 = ruEXT080.aru * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[4].plot(np.arange(2101,2200),tmp3,color=col080)
bias_2101_2120_4=np.mean(tmp3[0:20]-tmp2[86:106])
#-
tmp1 = ruEXT090.aru * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[4].plot(np.arange(2101,2200),tmp3,color=col090)
bias_2101_2120_5=np.mean(tmp3[0:20]-tmp2[86:106])
#-----
tmp1 = ruBAC050.aru * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[4].plot(np.arange(2015,2081),tmp3,color=col050)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[4].text(2025,3000,'r = 0.5  ('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120_1,precision=1)+')',fontsize=14,color=col050)
#-
tmp1 = ruBAC060.aru * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[4].plot(np.arange(2015,2081),tmp3,color=col060)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[4].text(2025,2600,'r = 0.6  ('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120_2,precision=1)+')',fontsize=14,color=col060)
#-
tmp1 = ruBAC070.aru * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[4].plot(np.arange(2015,2081),tmp3,color=col070)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[4].text(2025,2200,'r = 0.7  ('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120_3,precision=1)+')',fontsize=14,color=col070)
#-
tmp1 = ruBAC080.aru * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[4].plot(np.arange(2015,2081),tmp3,color=col080)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[4].text(2025,1800,'r = 0.8  ('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120_4,precision=1)+')',fontsize=14,color=col080)
#-
tmp1 = ruBAC090.aru * msk_grn * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[4].plot(np.arange(2015,2081),tmp3,color=col090)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[4].text(2025,1400,'r = 0.9  ('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120_5,precision=1)+')',fontsize=14,color=col090)
#-----
axs[4].set_title('(e) Anom. excess liquid water over the grounded i. s.',fontsize=16,weight='bold')
axs[4].tick_params(axis='both', labelsize=12)
axs[4].set_ylabel('Gt/yr',fontsize=14)
axs[4].set_xlim([2015,2200])
axs[4].set_ylim([-100,3400])

#----- f -----
axs[5].fill([2080.5,2100.5,2100.5,2080.5,2080.5],[-200,-200,4800,4800,-200],color='lightgrey',edgecolor=None)
axs[5].plot([1980,2200],[0,0],linewidth=0.5,color='k')
tmp1 = ruhis.aru * msk_isf * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[5].plot(np.arange(1980,2015),tmp2,color='black')
tmp1 = ru585.aru * msk_isf * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[5].plot(np.arange(2015,2200),tmp2,color='black')
#-
tmp1 = ruEXT050.aru * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[5].plot(np.arange(2101,2200),tmp3,color=col050)
bias_2101_2120_1=np.mean(tmp3[0:20]-tmp2[86:106])
#-
tmp1 = ruEXT060.aru * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[5].plot(np.arange(2101,2200),tmp3,color=col060)
bias_2101_2120_2=np.mean(tmp3[0:20]-tmp2[86:106])
#-
tmp1 = ruEXT070.aru * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[5].plot(np.arange(2101,2200),tmp3,color=col070)
bias_2101_2120_3=np.mean(tmp3[0:20]-tmp2[86:106])
#-
tmp1 = ruEXT080.aru * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[5].plot(np.arange(2101,2200),tmp3,color=col080)
bias_2101_2120_4=np.mean(tmp3[0:20]-tmp2[86:106])
#-
tmp1 = ruEXT090.aru * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,99)).rolling(min_periods=1,time=5,center=True).mean().values
axs[5].plot(np.arange(2101,2200),tmp3,color=col090)
bias_2101_2120_5=np.mean(tmp3[0:20]-tmp2[86:106])
#-----
tmp1 = ruBAC050.aru * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[5].plot(np.arange(2015,2081),tmp3,color=col050)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[5].text(2025,4000,'r = 0.5  ('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120_1,precision=1)+')',fontsize=14,color=col050)
#-
tmp1 = ruBAC060.aru * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[5].plot(np.arange(2015,2081),tmp3,color=col060)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[5].text(2025,3400,'r = 0.6  ('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120_2,precision=1)+')',fontsize=14,color=col060)
#-
tmp1 = ruBAC070.aru * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[5].plot(np.arange(2015,2081),tmp3,color=col070)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[5].text(2025,2800,'r = 0.7  ('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120_3,precision=1)+')',fontsize=14,color=col070)
#-
tmp1 = ruBAC080.aru * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[5].plot(np.arange(2015,2081),tmp3,color=col080)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[5].text(2025,2200,'r = 0.8  ('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120_4,precision=1)+')',fontsize=14,color=col080)
#-
tmp1 = ruBAC090.aru * msk_isf * fac; tmp3=tmp1.sum(dim=["x","y"]).isel(time=slice(0,66)).rolling(min_periods=1,time=5,center=True).mean().values
axs[5].plot(np.arange(2015,2081),tmp3,color=col090)
bias_2061_2080=np.mean(tmp3[46:66]-tmp2[46:66])
axs[5].text(2025,1600,'r = 0.9  ('+np.array2string(bias_2061_2080,precision=1)+' / '+np.array2string(bias_2101_2120_5,precision=1)+')',fontsize=14,color=col090)
#-----
axs[5].set_title('(f) Anom. excess liquid water over the ice shelves',fontsize=16,weight='bold')
axs[5].tick_params(axis='both', labelsize=12)
axs[5].set_ylabel('Gt/yr',fontsize=14)
axs[5].set_xlim([2015,2200])
axs[5].set_ylim([-200,4800])

figname='timeseries_extended_to_2200_IPSL_new.pdf'

fig.savefig(figname)
