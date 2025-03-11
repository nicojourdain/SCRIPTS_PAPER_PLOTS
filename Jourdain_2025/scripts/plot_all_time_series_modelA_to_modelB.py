import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

fig, axs = plt.subplots(nrows=3,ncols=2,figsize=(18.0,27.0))
axs = axs.ravel()

#modelB='CESM2' ; memberB='r11i1p1f1' ; scenarB='ssp585'
#modelB='UKESM1-0-LL' ; memberB='r1i1p1f2' ; scenarB='ssp585'
#modelB='CNRM-CM6-1' ; memberB='r1i1p1f2' ; scenarB='ssp585'
#modelB='MPI-ESM1-2-HR' ; memberB='r1i1p1f1' ; scenarB='ssp585'
modelB='ACCESS1-3' ; memberB='r1i1p1' ; scenarB='rcp85'

modelA  = [ 'UKESM1-0-LL', 'IPSL-CM6A-LR' , 'MPI-ESM1-2-HR' , 'CNRM-CM6-1', 'CESM2'     , 'ACCESS1-3' , 'NorESM1-M' ]
memberA = [ 'r1i1p1f2'   , 'r1i1p1f1'     , 'r1i1p1f1'      , 'r1i1p1f2'  , 'r11i1p1f1' , 'r1i1p1'    , 'r1i1p1'    ]
scenarA = [ 'ssp585'     , 'ssp585'       , 'ssp585'        , 'ssp585'    , 'ssp585'    , 'rcp85'     , 'rcp85'     ]
colorA  = [ 'darkgrey'   , 'fuchsia'      , 'mediumblue'    , 'orangered' , 'gold'      , 'olivedrab' , 'lightpink' ]
NmodA=np.size(modelA)

# input files
grd=xr.open_dataset('RCM_ice_regrid_04000m.nc2',decode_cf=False)
#-
smbhisB=xr.open_dataset('MAR-'+modelB+'_asmb_1980-2014_histo_regrid_04000m.nc',decode_cf=False)
mehisB=xr.open_dataset('MAR-'+modelB+'_ame_1980-2014_histo_regrid_04000m.nc',decode_cf=False)
ruhisB=xr.open_dataset('MAR-'+modelB+'_aru_1980-2014_histo_regrid_04000m.nc',decode_cf=False)
#-
smbfutB=xr.open_dataset('MAR-'+modelB+'_asmb_2015-2100_'+scenarB+'_regrid_04000m.nc',decode_cf=False)
mefutB=xr.open_dataset('MAR-'+modelB+'_ame_2015-2100_'+scenarB+'_regrid_04000m.nc',decode_cf=False)
rufutB=xr.open_dataset('MAR-'+modelB+'_aru_2015-2100_'+scenarB+'_regrid_04000m.nc',decode_cf=False)
#-
for kmodA in np.arange(NmodA):
  for var in ['smb', 'me', 'ru']:
    print('MAR-"+modelB+"-"+memberB+"_a"+var+"_1980-2014_histo_regrid_04000m_FROM_"+modelA[kmodA]+"-"+memberA[kmodA]+"-histo.nc')
    exec(var+'hisA'+kmodA.astype('str')+" = xr.open_dataset('MAR-"+modelB+"-"+memberB+"_a"+var+"_1980-2014_histo_regrid_04000m_FROM_"+modelA[kmodA]+"-"+memberA[kmodA]+"-histo.nc',decode_cf=False)")
    exec(var+"futA"+kmodA.astype('str')+" = xr.open_dataset('MAR-"+modelB+"-"+memberB+"_a"+var+"_2015-2100_"+scenarB+"_regrid_04000m_FROM_"+modelA[kmodA]+"-"+memberA[kmodA]+"-"+scenarA[kmodA]+".nc',decode_cf=False)")

# ice shelf and grounded ice sheet masks :
msk_isf = ( grd.ICE_MAR - grd.GROUND ) * grd.af2
msk_grn = grd.GROUND * grd.af2

# conversion to Gt/yr
fac=4.e3*4.e3*1.e-12*86400*365.25

bias_smb_grn_2015_2100=np.zeros((NmodA))
bias_me_grn_2015_2100=np.zeros((NmodA))
bias_ru_grn_2015_2100=np.zeros((NmodA))
bias_smb_isf_2015_2100=np.zeros((NmodA))
bias_me_isf_2015_2100=np.zeros((NmodA))
bias_ru_isf_2015_2100=np.zeros((NmodA))

#----- a -----
axs[0].plot([1980,2100],[0,0],linewidth=0.5,color='grey')
tmp1 = smbhisB.asmb * msk_grn * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[0].plot(np.arange(1980,2015),tmp2,color='black')
tmp1 = smbfutB.asmb * msk_grn * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[0].plot(np.arange(2015,2101),tmp2,color='black')
mean_smb_grn_2015_2100=np.mean(tmp2)
#-
for kmodA in np.arange(NmodA):
  tmp1 = eval("smbhisA"+kmodA.astype('str')).asmb * msk_grn * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
  axs[0].plot(np.arange(1980,2015),tmp2,color=colorA[kmodA])
  tmp1 = eval("smbfutA"+kmodA.astype('str')).asmb * msk_grn * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
  axs[0].plot(np.arange(2015,2101),tmp2,color=colorA[kmodA])
  bias_smb_grn_2015_2100[kmodA]=np.mean(tmp2)-mean_smb_grn_2015_2100
#-
ymin, ymax = axs[0].get_ylim()
dy=ymax-ymin
axs[0].text(1985,ymax-0.05*dy,'Original MAR-'+modelB+'  (mean='+np.array2string(mean_smb_grn_2015_2100,precision=1)+')',fontsize=12,color='black')
for kmodA in np.arange(NmodA):
  axs[0].text(1985,ymax-0.05*(kmodA+2)*dy,'MAR-'+modelB+' from MAR-'+modelA[kmodA]+' (bias='+np.array2string(bias_smb_grn_2015_2100[kmodA],precision=1)+')',fontsize=12,color=colorA[kmodA])
axs[0].set_title('(a) Grounded Ice Sheet SMB Anomaly',fontsize=16,fontweight='bold')
axs[0].tick_params(axis='both', labelsize=12)
axs[0].set_ylabel('Gt/yr',fontsize=14)
axs[0].set_xlim([1980,2100])

#----- b -----
axs[1].plot([1980,2100],[0,0],linewidth=0.5,color='grey')
tmp1 = smbhisB.asmb * msk_isf * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[1].plot(np.arange(1980,2015),tmp2,color='black')
tmp1 = smbfutB.asmb * msk_isf * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[1].plot(np.arange(2015,2101),tmp2,color='black')
mean_smb_isf_2015_2100=np.mean(tmp2)
#-
for kmodA in np.arange(NmodA):
  tmp1 = eval("smbhisA"+kmodA.astype('str')).asmb * msk_isf * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
  axs[1].plot(np.arange(1980,2015),tmp2,color=colorA[kmodA])
  tmp1 = eval("smbfutA"+kmodA.astype('str')).asmb * msk_isf * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
  axs[1].plot(np.arange(2015,2101),tmp2,color=colorA[kmodA])
  bias_smb_isf_2015_2100[kmodA]=np.mean(tmp2)-mean_smb_isf_2015_2100
#-
ymin, ymax = axs[1].get_ylim()
dy=ymax-ymin
axs[1].text(1985,ymax-0.5*dy,'Original MAR-'+modelB+'  (mean='+np.array2string(mean_smb_isf_2015_2100,precision=1)+')',fontsize=12,color='black')
for kmodA in np.arange(NmodA):
  axs[1].text(1985,ymax-0.5*dy-0.05*(kmodA+1)*dy,'MAR-'+modelB+' from MAR-'+modelA[kmodA]+' (bias='+np.array2string(bias_smb_isf_2015_2100[kmodA],precision=1)+')',fontsize=12,color=colorA[kmodA])
axs[1].set_title('(b) Ice Shelves SMB Anomaly',fontsize=16,fontweight='bold')
axs[1].tick_params(axis='both', labelsize=12)
axs[1].set_ylabel('Gt/yr',fontsize=14)
axs[1].set_xlim([1980,2100])

#----- c -----
axs[2].plot([1980,2100],[0,0],linewidth=0.5,color='grey')
tmp1 = mehisB.ame * msk_grn * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[2].plot(np.arange(1980,2015),tmp2,color='black')
tmp1 = mefutB.ame * msk_grn * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[2].plot(np.arange(2015,2101),tmp2,color='black')
mean_me_grn_2015_2100=np.mean(tmp2)
#-
for kmodA in np.arange(NmodA):
  tmp1 = eval("mehisA"+kmodA.astype('str')).ame * msk_grn * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
  axs[2].plot(np.arange(1980,2015),tmp2,color=colorA[kmodA])
  tmp1 = eval("mefutA"+kmodA.astype('str')).ame * msk_grn * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
  axs[2].plot(np.arange(2015,2101),tmp2,color=colorA[kmodA])
  bias_me_grn_2015_2100[kmodA]=np.mean(tmp2)-mean_me_grn_2015_2100
#-
ymin, ymax = axs[2].get_ylim()
dy=ymax-ymin
axs[2].text(1985,ymax-0.05*dy,'Original MAR-'+modelB+'  (mean='+np.array2string(mean_me_grn_2015_2100,precision=1)+')',fontsize=12,color='black')
for kmodA in np.arange(NmodA):
  axs[2].text(1985,ymax-0.05*(kmodA+2)*dy,'MAR-'+modelB+' from MAR-'+modelA[kmodA]+' (bias='+np.array2string(bias_me_grn_2015_2100[kmodA],precision=1)+')',fontsize=12,color=colorA[kmodA])
axs[2].set_title('(c) Grounded Ice Sheet Melt Anomaly',fontsize=16,fontweight='bold')
axs[2].tick_params(axis='both', labelsize=12)
axs[2].set_ylabel('Gt/yr',fontsize=14)
axs[2].set_xlim([1980,2100])

#----- d -----
axs[3].plot([1980,2100],[0,0],linewidth=0.5,color='grey')
tmp1 = mehisB.ame * msk_isf * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[3].plot(np.arange(1980,2015),tmp2,color='black')
tmp1 = mefutB.ame * msk_isf * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[3].plot(np.arange(2015,2101),tmp2,color='black')
mean_me_isf_2015_2100=np.mean(tmp2)
#-
for kmodA in np.arange(NmodA):
  tmp1 = eval("mehisA"+kmodA.astype('str')).ame * msk_isf * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
  axs[3].plot(np.arange(1980,2015),tmp2,color=colorA[kmodA])
  tmp1 = eval("mefutA"+kmodA.astype('str')).ame * msk_isf * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
  axs[3].plot(np.arange(2015,2101),tmp2,color=colorA[kmodA])
  bias_me_isf_2015_2100[kmodA]=np.mean(tmp2)-mean_me_isf_2015_2100
#-
ymin, ymax = axs[3].get_ylim()
dy=ymax-ymin
axs[3].text(1985,ymax-0.05*dy,'Original MAR-'+modelB+'  (mean='+np.array2string(mean_me_isf_2015_2100,precision=1)+')',fontsize=12,color='black')
for kmodA in np.arange(NmodA):
  axs[3].text(1985,ymax-0.05*(kmodA+2)*dy,'MAR-'+modelB+' from MAR-'+modelA[kmodA]+' (bias='+np.array2string(bias_me_isf_2015_2100[kmodA],precision=1)+')',fontsize=12,color=colorA[kmodA])
axs[3].set_title('(d) Ice Shelves Melt Anomaly',fontsize=16,fontweight='bold')
axs[3].tick_params(axis='both', labelsize=12)
axs[3].set_ylabel('Gt/yr',fontsize=14)
axs[3].set_xlim([1980,2100])

#----- e -----
axs[4].plot([1980,2100],[0,0],linewidth=0.5,color='grey')
tmp1 = ruhisB.aru * msk_grn * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[4].plot(np.arange(1980,2015),tmp2,color='black')
tmp1 = rufutB.aru * msk_grn * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[4].plot(np.arange(2015,2101),tmp2,color='black')
mean_ru_grn_2015_2100=np.mean(tmp2)
#-
for kmodA in np.arange(NmodA):
  tmp1 = eval("ruhisA"+kmodA.astype('str')).aru * msk_grn * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
  axs[4].plot(np.arange(1980,2015),tmp2,color=colorA[kmodA])
  tmp1 = eval("rufutA"+kmodA.astype('str')).aru * msk_grn * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
  axs[4].plot(np.arange(2015,2101),tmp2,color=colorA[kmodA])
  bias_ru_grn_2015_2100[kmodA]=np.mean(tmp2)-mean_ru_grn_2015_2100
#-
ymin, ymax = axs[4].get_ylim()
dy=ymax-ymin
axs[4].text(1985,ymax-0.05*dy,'Original MAR-'+modelB+'  (mean='+np.array2string(mean_ru_grn_2015_2100,precision=1)+')',fontsize=12,color='black')
for kmodA in np.arange(NmodA):
  axs[4].text(1985,ymax-0.05*(kmodA+2)*dy,'MAR-'+modelB+' from MAR-'+modelA[kmodA]+' (bias='+np.array2string(bias_ru_grn_2015_2100[kmodA],precision=1)+')',fontsize=12,color=colorA[kmodA])
axs[4].set_title('(e) Grounded Ice Sheet Runoff Anomaly',fontsize=16,fontweight='bold')
axs[4].tick_params(axis='both', labelsize=12)
axs[4].set_ylabel('Gt/yr',fontsize=14)
axs[4].set_xlim([1980,2100])

#----- f -----
axs[5].plot([1980,2100],[0,0],linewidth=0.5,color='grey')
tmp1 = ruhisB.aru * msk_isf * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[5].plot(np.arange(1980,2015),tmp2,color='black')
tmp1 = rufutB.aru * msk_isf * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
axs[5].plot(np.arange(2015,2101),tmp2,color='black')
mean_ru_isf_2015_2100=np.mean(tmp2)
#-
for kmodA in np.arange(NmodA):
  tmp1 = eval("ruhisA"+kmodA.astype('str')).aru * msk_isf * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
  axs[5].plot(np.arange(1980,2015),tmp2,color=colorA[kmodA])
  tmp1 = eval("rufutA"+kmodA.astype('str')).aru * msk_isf * fac; tmp2=tmp1.sum(dim=["x","y"]).rolling(min_periods=1,time=5,center=True).mean().values
  axs[5].plot(np.arange(2015,2101),tmp2,color=colorA[kmodA])
  bias_ru_isf_2015_2100[kmodA]=np.mean(tmp2)-mean_ru_isf_2015_2100
#-
ymin, ymax = axs[5].get_ylim()
dy=ymax-ymin
axs[5].text(1985,ymax-0.05*dy,'Original MAR-'+modelB+'  (mean='+np.array2string(mean_ru_isf_2015_2100,precision=1)+')',fontsize=12,color='black')
for kmodA in np.arange(NmodA):
  axs[5].text(1985,ymax-0.05*(kmodA+2)*dy,'MAR-'+modelB+' from MAR-'+modelA[kmodA]+' (bias='+np.array2string(bias_ru_isf_2015_2100[kmodA],precision=1)+')',fontsize=12,color=colorA[kmodA])
axs[5].set_title('(f) Ice Shelves Runoff Anomaly',fontsize=16,fontweight='bold')
axs[5].tick_params(axis='both', labelsize=12)
axs[5].set_ylabel('Gt/yr',fontsize=14)
axs[5].set_xlim([1980,2100])

#-----

fileout='stats_other_models_to_'+modelB+'.npz'

np.savez(fileout,\
modelA=modelA,modelB=modelB,memberA=memberA,scenarA=scenarA,scenarB=scenarB,colorA=colorA,\
mean_smb_grn_2015_2100=mean_smb_grn_2015_2100,\
mean_smb_isf_2015_2100=mean_smb_isf_2015_2100,\
mean_me_grn_2015_2100=mean_me_grn_2015_2100,\
mean_me_isf_2015_2100=mean_me_isf_2015_2100,\
mean_ru_grn_2015_2100=mean_ru_grn_2015_2100,\
mean_ru_isf_2015_2100=mean_ru_isf_2015_2100,\
bias_smb_grn_2015_2100=bias_smb_grn_2015_2100,\
bias_smb_isf_2015_2100=bias_smb_isf_2015_2100,\
bias_me_grn_2015_2100=bias_me_grn_2015_2100,\
bias_me_isf_2015_2100=bias_me_isf_2015_2100,\
bias_ru_grn_2015_2100=bias_ru_grn_2015_2100,\
bias_ru_isf_2015_2100=bias_ru_isf_2015_2100)

figname='timeseries_other_models_to_'+modelB+'.pdf'

fig.savefig(figname)
