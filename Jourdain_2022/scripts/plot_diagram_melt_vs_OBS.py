import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors

#------------------------------------
# Define new colors:

col1='cornflowerblue'
colors.ColorConverter.colors['blue1'] = (colors.to_rgb('cornflowerblue')[0]*0.60,colors.to_rgb('cornflowerblue')[1]*0.60,colors.to_rgb('cornflowerblue')[2]*0.60)
col3='blue1'
colors.ColorConverter.colors['blue2'] = (colors.to_rgb('cornflowerblue')[0]*0.30,colors.to_rgb('cornflowerblue')[1]*0.30,colors.to_rgb('cornflowerblue')[2]*0.30)
col5='blue2'

#------------------------------------
# Observation-based estimates

## Naughten et al. (2021), Dutrieux et al. (2014)
#melt_PIG_obs=np.zeros((3,5))
#melt_PIG_obs[0,:] = np.array([ 41.7, 65.7, 59.6, 30.0, 22.2 ])
#melt_PIG_obs[1,:] = np.array([ 46.3, 72.9, 66.2, 33.2, 36.9 ])
#melt_PIG_obs[2,:] = np.array([ 50.9, 80.2, 72.8, 36.5, 51.6 ])
#year_PIG_obs = np.array([ 1994, 2009, 2010, 2012, 2014 ])
# Joughin et al. (2021):
melt_PIG_obs=np.zeros((3,6))
melt_PIG_obs[0,:] = np.array([ 37.7,  68.2, 57.2, 52.6, 11.7,  1.2 ])
melt_PIG_obs[1,:] = np.array([ 47.8,  99.0, 78.2, 60.4, 31.5, 30.6 ])
melt_PIG_obs[2,:] = np.array([ 57.7, 129.5, 99.0, 68.2, 51.4, 59.9 ])
year_PIG_obs = np.array([      1994,  2007, 2009, 2010, 2012, 2014 ])

# Jenkins et al. (2018)
melt_DOT_obs=np.zeros((3,8))
melt_DOT_obs[0,:] = np.array([ 17.0, 40.0,  9.0, 60.0, 42.0, 14.0, 13.0,  4.0 ])
melt_DOT_obs[1,:] = np.array([ 25.0, 55.0, 45.0, 90.0, 52.0, 20.0, 21.0, 20.0 ])
melt_DOT_obs[2,:] = np.array([ 35.0, 70.0, 80.0,123.0, 65.0, 26.0, 29.0, 34.0 ])
year_DOT_obs = np.array([      2000, 2006, 2007, 2009, 2011, 2012, 2014, 2016 ])


N=1000000

tmp_PIG=np.zeros((N))
for kk in np.arange(np.size(year_PIG_obs)):
  tmp_PIG = tmp_PIG + np.random.normal(loc=melt_PIG_obs[1,kk], scale=0.5*(melt_PIG_obs[2,kk]-melt_PIG_obs[0,kk]), size=N)
tmp_PIG = tmp_PIG / np.size(year_PIG_obs)
tmp_PIG = tmp_PIG*1.e9 / 6000.e6 # Gt/yr -> m.w.e./yr
mean_PIG_oce = np.mean(tmp_PIG)
std_PIG_oce = np.std(tmp_PIG)
p975_PIG_oce = np.percentile(tmp_PIG,97.5)
p025_PIG_oce = np.percentile(tmp_PIG,2.5)

tmp_DOT=np.zeros((N))
for kk in np.arange(np.size(year_DOT_obs)):
  tmp_DOT = tmp_DOT + np.random.normal(loc=melt_DOT_obs[1,kk], scale=0.5*(melt_DOT_obs[2,kk]-melt_DOT_obs[0,kk]), size=N)
print(np.shape(tmp_DOT))
tmp_DOT = tmp_DOT / np.size(year_DOT_obs)
tmp_DOT = tmp_DOT*1.e9 / 5700.e6 # Gt/yr -> m.w.e./yr
mean_DOT_oce = np.mean(tmp_DOT)
std_DOT_oce = np.std(tmp_DOT)
p975_DOT_oce = np.percentile(tmp_DOT,97.5)
p025_DOT_oce = np.percentile(tmp_DOT,2.5)

# Adusumili (2019), for 1994-2018, converted from meters of ice to m.w.e. :
mean_VEN_Adu =  5.1 * 0.920 ; p95_VEN_Adu =  2.0 * 0.920 # Venable
mean_ABB_Adu =  1.5 * 0.920 ; p95_ABB_Adu =  1.5 * 0.920 # Abbot
mean_COS_Adu =  1.0 * 0.920 ; p95_COS_Adu =  1.5 * 0.920 # Cosgrove
mean_PIG_Adu = 14.0 * 0.920 ; p95_PIG_Adu =  1.6 * 0.920 # PIG
mean_THW_Adu = 26.7 * 0.920 ; p95_THW_Adu =  2.4 * 0.920 # Thwaites
mean_CRO_Adu =  7.8 * 0.920 ; p95_CRO_Adu =  1.8 * 0.920 # Crosson
mean_DOT_Adu =  5.4 * 0.920 ; p95_DOT_Adu =  1.6 * 0.920 # Dotson
mean_GET_Adu =  4.2 * 0.920 ; p95_GET_Adu =  1.4 * 0.920 # Getz

# Rignot et al. (2013), for ~2003-2008, converting from ±std to 95% interval
mean_VEN_Rig =  6.1 ; p95_VEN_Rig =  0.7 * 2.0 # Venable
mean_ABB_Rig =  1.7 ; p95_ABB_Rig =  0.6 * 2.0 # Abbot
mean_COS_Rig =  2.8 ; p95_COS_Rig =  0.7 * 2.0 # Cosgrove
mean_PIG_Rig = 16.2 ; p95_PIG_Rig =  1.0 * 2.0 # PIG
mean_THW_Rig = 17.7 ; p95_THW_Rig =  1.0 * 2.0 # Thwaites
mean_CRO_Rig = 11.9 ; p95_CRO_Rig =  1.0 * 2.0 # Crosson
mean_DOT_Rig =  7.8 ; p95_DOT_Rig =  0.6 * 2.0 # Dotson
mean_GET_Rig =  4.3 ; p95_GET_Rig =  0.4 * 2.0 # Getz

# Depoorter et al. (2013), for ~2003-2008, converting from ±std to 95% interval
mean_VEN_Dep =  4.82; p95_VEN_Dep =  0.95* 2.0 # Venable
mean_ABB_Dep =  2.72; p95_ABB_Dep =  0.70* 2.0 # Abbot
mean_COS_Dep =  3.74; p95_COS_Dep =  0.89* 2.0 # Cosgrove
mean_PIG_Dep = 15.96; p95_PIG_Dep =  2.38* 2.0 # PIG
mean_THW_Dep = 15.22; p95_THW_Dep =  3.87* 2.0 # Thwaites
mean_GET_Dep =  4.09; p95_GET_Dep =  0.68* 2.0 # Getz

#---------------------------------------------------------------------------------

npzfile1 = np.load('../OUTPUT_AMU/melt_AMUXL12_GNJ002_BM02MAR.npz')
npzfile3 = np.load('../OUTPUT_AMU/melt_AMUXL12_GNJ002_BM03MAR.npz')
npzfile5 = np.load('../OUTPUT_AMU/melt_AMUXL12_GNJ002_BM04MAR.npz')

time1=npzfile1['time']; mt1=np.size(time1)
time3=npzfile3['time']; mt3=np.size(time3)
time5=npzfile5['time']; mt5=np.size(time5)

melt_mw_yr1=npzfile1['melt_mw_yr']
melt_mw_yr3=npzfile3['melt_mw_yr']
melt_mw_yr5=npzfile5['melt_mw_yr']

name_isf1=npzfile1['name_isf']
name_isf3=npzfile3['name_isf']
name_isf5=npzfile5['name_isf']

for kk in np.arange(0,np.size(name_isf3),1):
  if "Venable" in name_isf3[kk]:
    k_VEN=kk
  elif "Abbot" in name_isf3[kk]:
    k_ABB=kk
  elif "Cosgrove" in name_isf3[kk]:
    k_COS=kk
  elif "Pine Island" in name_isf3[kk]:
    k_PIG=kk
  elif "Thwaites" in name_isf3[kk]:
    k_THW=kk
  elif "Crosson" in name_isf3[kk]:
    k_CRO=kk
  elif "Dotson" in name_isf3[kk]:
    k_DOT=kk
  elif "Getz" in name_isf3[kk]:
    k_GET=kk

#-----------------------------------------------
# PLOT:

def xrect(xmid,width):
   x1=xmid-width*0.5
   x2=xmid+width*0.5
   x=[x1,x2,x2,x1,x1]
   return x

def yrect(ymax):
   y=[0.0,0.0,ymax,ymax,0.0]
   return y

def ebar(x,y,low,high,thick=False,col='grey'):
   dbar=0.1
   if thick:
     lwdth=1.2
     ma='s'
   else:
     lwdth=0.5
     ma='d'
   plt.plot([x, x],[low, high],color=col,linewidth=lwdth,zorder=1)
   plt.plot([x-0.5*dbar, x+0.5*dbar],[low, low],color=col,linewidth=lwdth,zorder=2)
   plt.plot([x-0.5*dbar, x+0.5*dbar],[high, high],color=col,linewidth=lwdth,zorder=3)
   plt.scatter(x,y,s=5,c=col,marker=ma,zorder=4)

fig, ax = plt.subplots()

dplt=0.90

# Getz:
xplt=1
xti=np.array([xplt])
plt.fill(xrect(xplt-0.2,0.20),yrect(np.mean(melt_mw_yr1[k_GET,120:361])),col1)
ebar(xplt-0.2,np.mean(melt_mw_yr1[k_GET,120:361]),np.percentile(melt_mw_yr1[k_GET,120:361],2.5),np.percentile(melt_mw_yr1[k_GET,120:361],97.5))
plt.fill(xrect(xplt+0.0,0.20),yrect(np.mean(melt_mw_yr3[k_GET,120:361])),col3)
ebar(xplt+0.0,np.mean(melt_mw_yr3[k_GET,120:361]),np.percentile(melt_mw_yr3[k_GET,120:361],2.5),np.percentile(melt_mw_yr3[k_GET,120:361],97.5))
plt.fill(xrect(xplt+0.2,0.20),yrect(np.mean(melt_mw_yr5[k_GET,120:361])),col5)
ebar(xplt+0.2,np.mean(melt_mw_yr5[k_GET,120:361]),np.percentile(melt_mw_yr5[k_GET,120:361],2.5),np.percentile(melt_mw_yr5[k_GET,120:361],97.5))
ebar(xplt-0.03,mean_GET_Adu,mean_GET_Adu-p95_GET_Adu,mean_GET_Adu+p95_GET_Adu,thick=True,col='red')
ebar(xplt+0.03,mean_GET_Rig,mean_GET_Rig-p95_GET_Rig,mean_GET_Rig+p95_GET_Rig,thick=True,col='purple')
ebar(xplt+0.09,mean_GET_Dep,mean_GET_Dep-p95_GET_Dep,mean_GET_Dep+p95_GET_Dep,thick=True,col='cyan')

# Dotson:
xplt=xplt+dplt
xti=np.append(xti,np.array([xplt]))
plt.fill(xrect(xplt-0.2,0.20),yrect(np.mean(melt_mw_yr1[k_DOT,120:361])),col1)
ebar(xplt-0.2,np.mean(melt_mw_yr1[k_DOT,120:361]),np.percentile(melt_mw_yr1[k_DOT,120:361],2.5),np.percentile(melt_mw_yr1[k_DOT,120:361],97.5))
plt.fill(xrect(xplt+0.0,0.20),yrect(np.mean(melt_mw_yr3[k_DOT,120:361])),col3)
ebar(xplt+0.0,np.mean(melt_mw_yr3[k_DOT,120:361]),np.percentile(melt_mw_yr3[k_DOT,120:361],2.5),np.percentile(melt_mw_yr3[k_DOT,120:361],97.5))
plt.fill(xrect(xplt+0.2,0.20),yrect(np.mean(melt_mw_yr5[k_DOT,120:361])),col5)
ebar(xplt+0.2,np.mean(melt_mw_yr5[k_DOT,120:361]),np.percentile(melt_mw_yr5[k_DOT,120:361],2.5),np.percentile(melt_mw_yr5[k_DOT,120:361],97.5))
ebar(xplt-0.09,mean_DOT_oce,p025_DOT_oce,p975_DOT_oce,thick=True,col='y')
ebar(xplt-0.03,mean_DOT_Adu,mean_DOT_Adu-p95_DOT_Adu,mean_DOT_Adu+p95_DOT_Adu,thick=True,col='red')
ebar(xplt+0.03,mean_DOT_Rig,mean_DOT_Rig-p95_DOT_Rig,mean_DOT_Rig+p95_DOT_Rig,thick=True,col='purple')

# Crosson:
xplt=xplt+dplt
xti=np.append(xti,np.array([xplt]))
plt.fill(xrect(xplt-0.2,0.20),yrect(np.mean(melt_mw_yr1[k_CRO,120:361])),col1)
ebar(xplt-0.2,np.mean(melt_mw_yr1[k_CRO,120:361]),np.percentile(melt_mw_yr1[k_CRO,120:361],2.5),np.percentile(melt_mw_yr1[k_CRO,120:361],97.5))
plt.fill(xrect(xplt+0.0,0.20),yrect(np.mean(melt_mw_yr3[k_CRO,120:361])),col3)
ebar(xplt+0.0,np.mean(melt_mw_yr3[k_CRO,120:361]),np.percentile(melt_mw_yr3[k_CRO,120:361],2.5),np.percentile(melt_mw_yr3[k_CRO,120:361],97.5))
plt.fill(xrect(xplt+0.2,0.20),yrect(np.mean(melt_mw_yr5[k_CRO,120:361])),col5)
ebar(xplt+0.2,np.mean(melt_mw_yr5[k_CRO,120:361]),np.percentile(melt_mw_yr5[k_CRO,120:361],2.5),np.percentile(melt_mw_yr5[k_CRO,120:361],97.5))
ebar(xplt-0.03,mean_CRO_Adu,mean_CRO_Adu-p95_CRO_Adu,mean_CRO_Adu+p95_CRO_Adu,thick=True,col='red')
ebar(xplt+0.03,mean_CRO_Rig,mean_CRO_Rig-p95_CRO_Rig,mean_CRO_Rig+p95_CRO_Rig,thick=True,col='purple')

# Thwaites:
xplt=xplt+dplt
xti=np.append(xti,np.array([xplt]))
plt.fill(xrect(xplt-0.2,0.20),yrect(np.mean(melt_mw_yr1[k_THW,120:361])),col1)
ebar(xplt-0.2,np.mean(melt_mw_yr1[k_THW,120:361]),np.percentile(melt_mw_yr1[k_THW,120:361],2.5),np.percentile(melt_mw_yr1[k_THW,120:361],97.5))
plt.fill(xrect(xplt+0.0,0.20),yrect(np.mean(melt_mw_yr3[k_THW,120:361])),col3)
ebar(xplt+0.0,np.mean(melt_mw_yr3[k_THW,120:361]),np.percentile(melt_mw_yr3[k_THW,120:361],2.5),np.percentile(melt_mw_yr3[k_THW,120:361],97.5))
plt.fill(xrect(xplt+0.2,0.20),yrect(np.mean(melt_mw_yr5[k_THW,120:361])),col5)
ebar(xplt+0.2,np.mean(melt_mw_yr5[k_THW,120:361]),np.percentile(melt_mw_yr5[k_THW,120:361],2.5),np.percentile(melt_mw_yr5[k_THW,120:361],97.5))
ebar(xplt-0.03,mean_THW_Adu,mean_THW_Adu-p95_THW_Adu,mean_THW_Adu+p95_THW_Adu,thick=True,col='red')
ebar(xplt+0.03,mean_THW_Rig,mean_THW_Rig-p95_THW_Rig,mean_THW_Rig+p95_THW_Rig,thick=True,col='purple')
ebar(xplt+0.09,mean_THW_Dep,mean_THW_Dep-p95_THW_Dep,mean_THW_Dep+p95_THW_Dep,thick=True,col='cyan')

# Pine Island:
xplt=xplt+dplt
xti=np.append(xti,np.array([xplt]))
plt.fill(xrect(xplt-0.2,0.20),yrect(np.mean(melt_mw_yr1[k_PIG,120:361])),col1)
ebar(xplt-0.2,np.mean(melt_mw_yr1[k_PIG,120:361]),np.percentile(melt_mw_yr1[k_PIG,120:361],2.5),np.percentile(melt_mw_yr1[k_PIG,120:361],97.5))
plt.fill(xrect(xplt+0.0,0.20),yrect(np.mean(melt_mw_yr3[k_PIG,120:361])),col3)
ebar(xplt+0.0,np.mean(melt_mw_yr3[k_PIG,120:361]),np.percentile(melt_mw_yr3[k_PIG,120:361],2.5),np.percentile(melt_mw_yr3[k_PIG,120:361],97.5))
plt.fill(xrect(xplt+0.2,0.20),yrect(np.mean(melt_mw_yr5[k_PIG,120:361])),col5)
ebar(xplt+0.2,np.mean(melt_mw_yr5[k_PIG,120:361]),np.percentile(melt_mw_yr5[k_PIG,120:361],2.5),np.percentile(melt_mw_yr5[k_PIG,120:361],97.5))
ebar(xplt-0.09,mean_PIG_oce,p025_PIG_oce,p975_PIG_oce,thick=True,col='y')
ebar(xplt-0.03,mean_PIG_Adu,mean_PIG_Adu-p95_PIG_Adu,mean_PIG_Adu+p95_PIG_Adu,thick=True,col='red')
ebar(xplt+0.03,mean_PIG_Rig,mean_PIG_Rig-p95_PIG_Rig,mean_PIG_Rig+p95_PIG_Rig,thick=True,col='purple')
ebar(xplt+0.09,mean_PIG_Dep,mean_PIG_Dep-p95_PIG_Dep,mean_PIG_Dep+p95_PIG_Dep,thick=True,col='cyan')

# Cosgrove:
xplt=xplt+dplt
xti=np.append(xti,np.array([xplt]))
plt.fill(xrect(xplt-0.2,0.20),yrect(np.mean(melt_mw_yr1[k_COS,120:361])),col1)
ebar(xplt-0.2,np.mean(melt_mw_yr1[k_COS,120:361]),np.percentile(melt_mw_yr1[k_COS,120:361],2.5),np.percentile(melt_mw_yr1[k_COS,120:361],97.5))
plt.fill(xrect(xplt+0.0,0.20),yrect(np.mean(melt_mw_yr3[k_COS,120:361])),col3)
ebar(xplt+0.0,np.mean(melt_mw_yr3[k_COS,120:361]),np.percentile(melt_mw_yr3[k_COS,120:361],2.5),np.percentile(melt_mw_yr3[k_COS,120:361],97.5))
plt.fill(xrect(xplt+0.2,0.20),yrect(np.mean(melt_mw_yr5[k_COS,120:361])),col5)
ebar(xplt+0.2,np.mean(melt_mw_yr5[k_COS,120:361]),np.percentile(melt_mw_yr5[k_COS,120:361],2.5),np.percentile(melt_mw_yr5[k_COS,120:361],97.5))
ebar(xplt-0.03,mean_COS_Adu,mean_COS_Adu-p95_COS_Adu,mean_COS_Adu+p95_COS_Adu,thick=True,col='red')
ebar(xplt+0.03,mean_COS_Rig,mean_COS_Rig-p95_COS_Rig,mean_COS_Rig+p95_COS_Rig,thick=True,col='purple')
ebar(xplt+0.09,mean_COS_Dep,mean_COS_Dep-p95_COS_Dep,mean_COS_Dep+p95_COS_Dep,thick=True,col='cyan')

# Abbot:
xplt=xplt+dplt
xti=np.append(xti,np.array([xplt]))
plt.fill(xrect(xplt-0.2,0.20),yrect(np.mean(melt_mw_yr1[k_ABB,120:361])),col1)
ebar(xplt-0.2,np.mean(melt_mw_yr1[k_ABB,120:361]),np.percentile(melt_mw_yr1[k_ABB,120:361],2.5),np.percentile(melt_mw_yr1[k_ABB,120:361],97.5))
plt.fill(xrect(xplt+0.0,0.20),yrect(np.mean(melt_mw_yr3[k_ABB,120:361])),col3)
ebar(xplt+0.0,np.mean(melt_mw_yr3[k_ABB,120:361]),np.percentile(melt_mw_yr3[k_ABB,120:361],2.5),np.percentile(melt_mw_yr3[k_ABB,120:361],97.5))
plt.fill(xrect(xplt+0.2,0.20),yrect(np.mean(melt_mw_yr5[k_ABB,120:361])),col5)
ebar(xplt+0.2,np.mean(melt_mw_yr5[k_ABB,120:361]),np.percentile(melt_mw_yr5[k_ABB,120:361],2.5),np.percentile(melt_mw_yr5[k_ABB,120:361],97.5))
ebar(xplt-0.03,mean_ABB_Adu,mean_ABB_Adu-p95_ABB_Adu,mean_ABB_Adu+p95_ABB_Adu,thick=True,col='red')
ebar(xplt+0.03,mean_ABB_Rig,mean_ABB_Rig-p95_ABB_Rig,mean_ABB_Rig+p95_ABB_Rig,thick=True,col='purple')
ebar(xplt+0.09,mean_ABB_Dep,mean_ABB_Dep-p95_ABB_Dep,mean_ABB_Dep+p95_ABB_Dep,thick=True,col='cyan')

# Venable:
xplt=xplt+dplt
xti=np.append(xti,np.array([xplt]))
plt.fill(xrect(xplt-0.2,0.20),yrect(np.mean(melt_mw_yr1[k_VEN,120:361])),col1)
ebar(xplt-0.2,np.mean(melt_mw_yr1[k_VEN,120:361]),np.percentile(melt_mw_yr1[k_VEN,120:361],2.5),np.percentile(melt_mw_yr1[k_VEN,120:361],97.5))
plt.fill(xrect(xplt+0.0,0.20),yrect(np.mean(melt_mw_yr3[k_VEN,120:361])),col3)
ebar(xplt+0.0,np.mean(melt_mw_yr3[k_VEN,120:361]),np.percentile(melt_mw_yr3[k_VEN,120:361],2.5),np.percentile(melt_mw_yr3[k_VEN,120:361],97.5))
plt.fill(xrect(xplt+0.2,0.20),yrect(np.mean(melt_mw_yr5[k_VEN,120:361])),col5)
ebar(xplt+0.2,np.mean(melt_mw_yr5[k_VEN,120:361]),np.percentile(melt_mw_yr5[k_VEN,120:361],2.5),np.percentile(melt_mw_yr5[k_VEN,120:361],97.5))
ebar(xplt-0.03,mean_VEN_Adu,mean_VEN_Adu-p95_VEN_Adu,mean_VEN_Adu+p95_VEN_Adu,thick=True,col='red')
ebar(xplt+0.03,mean_VEN_Rig,mean_VEN_Rig-p95_VEN_Rig,mean_VEN_Rig+p95_VEN_Rig,thick=True,col='purple')
ebar(xplt+0.09,mean_VEN_Dep,mean_VEN_Dep-p95_VEN_Dep,mean_VEN_Dep+p95_VEN_Dep,thick=True,col='cyan')

# Legend:
sq=[0.,0.,1.,1.,0.]
plt.fill(xrect(0.7*xplt,0.3),33.*np.ones((5))+sq,col1)
plt.text(0.7*xplt+0.35,33.5,'NEMO 1989-2009 (A)',fontsize=8,HorizontalAlignment='left',VerticalAlignment='center')
plt.fill(xrect(0.7*xplt,0.3),29.*np.ones((5))+sq,col3)
plt.text(0.7*xplt+0.35,29.5,'NEMO 1989-2009 (B)',fontsize=8,HorizontalAlignment='left',VerticalAlignment='center')
plt.fill(xrect(0.7*xplt,0.3),25.*np.ones((5))+sq,col5)
plt.text(0.7*xplt+0.35,25.5,'NEMO 1989-2009 (C)',fontsize=8,HorizontalAlignment='left',VerticalAlignment='center')
ebar(0.1*xplt,33.5,33,34,thick=True,col='y')
plt.text(0.13*xplt,33.5,'oceanographic estimates',fontsize=8,HorizontalAlignment='left',VerticalAlignment='center')
ebar(0.1*xplt,31.5,31,32,thick=True,col='red')
plt.text(0.13*xplt,31.5,'Adusumilli et al. (2020)',fontsize=8,HorizontalAlignment='left',VerticalAlignment='center')
ebar(0.1*xplt,29.5,29,30,thick=True,col='purple')
plt.text(0.13*xplt,29.5,'Rignot et al. (2013)',fontsize=8,HorizontalAlignment='left',VerticalAlignment='center')
ebar(0.1*xplt,27.5,27,28,thick=True,col='cyan')
plt.text(0.13*xplt,27.5,'Depoorter et al. (2013)',fontsize=8,HorizontalAlignment='left',VerticalAlignment='center')

plt.ylim(bottom=0.0)
plt.ylabel('ice-shelf melt rate (m.w.e./yr)',size=8)
plt.yticks(fontsize=8)
ax.set_xticks(xti)
ax.set_xticklabels(('Getz','Dotson','Crosson','Thwaites','Pine Isl.','Cosgrove','Abbot','Venable'),size=8)

#plt.show()

fig.savefig('diagram_melt_vs_OBS.jpg')
fig.savefig('diagram_melt_vs_OBS.pdf')
