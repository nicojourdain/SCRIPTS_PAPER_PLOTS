import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors

#------------------------------------
# Define new colors:

col1='cornflowerblue'
col2='orange'
colors.ColorConverter.colors['blue1'] = (colors.to_rgb('cornflowerblue')[0]*0.60,colors.to_rgb('cornflowerblue')[1]*0.60,colors.to_rgb('cornflowerblue')[2]*0.60)
colors.ColorConverter.colors['orange1'] = (colors.to_rgb('orange')[0]*0.60,colors.to_rgb('orange')[1]*0.60,colors.to_rgb('orange')[2]*0.60)
col3='blue1'
col4='orange1'
colors.ColorConverter.colors['blue2'] = (colors.to_rgb('cornflowerblue')[0]*0.30,colors.to_rgb('cornflowerblue')[1]*0.30,colors.to_rgb('cornflowerblue')[2]*0.30)
colors.ColorConverter.colors['orange2'] = (colors.to_rgb('orange')[0]*0.30,colors.to_rgb('orange')[1]*0.30,colors.to_rgb('orange')[2]*0.30)
col5='blue2'
col6='orange2'

#---------------------------------------------------------------------------------

npzfile1 = np.load('../OUTPUT_AMU/melt_AMUXL12_GNJ002_BM02MAR.npz')
npzfile2 = np.load('../OUTPUT_AMU/melt_AMUXL12_GNJ002_BM02MARrcp85.npz')
npzfile3 = np.load('../OUTPUT_AMU/melt_AMUXL12_GNJ002_BM03MAR.npz')
npzfile4 = np.load('../OUTPUT_AMU/melt_AMUXL12_GNJ002_BM03MARrcBDY.npz')
npzfile5 = np.load('../OUTPUT_AMU/melt_AMUXL12_GNJ002_BM04MAR.npz')
npzfile6 = np.load('../OUTPUT_AMU/melt_AMUXL12_GNJ002_BM04MARrcp85.npz')
npzfile7 = np.load('../OUTPUT_AMU/melt_AMUXL12_GNJ002_BM03MARrcp85.npz') # no BDY perturbation


melt_mw_yr1=npzfile1['melt_mw_yr']
melt_mw_yr2=npzfile2['melt_mw_yr']
melt_mw_yr3=npzfile3['melt_mw_yr']
melt_mw_yr4=npzfile4['melt_mw_yr']
melt_mw_yr5=npzfile5['melt_mw_yr']
melt_mw_yr6=npzfile6['melt_mw_yr']
melt_mw_yr7=npzfile7['melt_mw_yr']

name_isf1=npzfile1['name_isf']
name_isf2=npzfile2['name_isf']
name_isf3=npzfile3['name_isf']
name_isf4=npzfile4['name_isf']
name_isf5=npzfile5['name_isf']
name_isf6=npzfile6['name_isf']
name_isf7=npzfile7['name_isf']

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

def xtri(xmid,width):
   x1=xmid
   x2=xmid+width
   x=[x1,x1,x2,x1]
   return x

def ytri(ymax,width):
   y=[ymax-0.5*width, ymax+0.5*width, ymax, ymax-0.5*width]
   return y

def ebar(x,y,low,high):
   dbar=0.1
   plt.plot([x, x],[low, high],color='grey',linewidth=0.5,zorder=91)
   plt.plot([x-0.5*dbar, x+0.5*dbar],[low, low],color='grey',linewidth=0.5,zorder=92)
   plt.plot([x-0.5*dbar, x+0.5*dbar],[high, high],color='grey',linewidth=0.5,zorder=93)
   plt.scatter(x,y,s=5,c='grey',marker='d',zorder=94)

fig, ax = plt.subplots()

dplt=1.2

#----------
# Getz:
xplt=1
xti=np.array([xplt])

plt.fill(xrect(xplt-0.25,0.2),yrect(np.mean(melt_mw_yr2[k_GET,120:361])),col2,zorder=11)
ebar(xplt-0.25,np.mean(melt_mw_yr2[k_GET,120:361]),np.percentile(melt_mw_yr2[k_GET,120:361],2.5),np.percentile(melt_mw_yr2[k_GET,120:361],97.5))
plt.fill(xrect(xplt-0.30,0.2),yrect(np.mean(melt_mw_yr1[k_GET,120:361])),col1,zorder=12)
ebar(xplt-0.30,np.mean(melt_mw_yr1[k_GET,120:361]),np.percentile(melt_mw_yr1[k_GET,120:361],2.5),np.percentile(melt_mw_yr1[k_GET,120:361],97.5))

plt.fill(xrect(xplt+0.025,0.2),yrect(np.mean(melt_mw_yr4[k_GET,120:361])),col4,zorder=13)
plt.plot([xplt-0.076,xplt+0.126],[np.mean(melt_mw_yr7[k_GET,120:361]),np.mean(melt_mw_yr7[k_GET,120:361])],color='white',linewidth=2.5,zorder=14)
ebar(xplt+0.025,np.mean(melt_mw_yr4[k_GET,120:361]),np.percentile(melt_mw_yr4[k_GET,120:361],2.5),np.percentile(melt_mw_yr4[k_GET,120:361],97.5))
plt.fill(xrect(xplt-0.025,0.2),yrect(np.mean(melt_mw_yr3[k_GET,120:361])),col3,zorder=15)
ebar(xplt-0.025,np.mean(melt_mw_yr3[k_GET,120:361]),np.percentile(melt_mw_yr3[k_GET,120:361],2.5),np.percentile(melt_mw_yr3[k_GET,120:361],97.5))

plt.fill(xrect(xplt+0.30,0.2),yrect(np.mean(melt_mw_yr6[k_GET,120:361])),col6,zorder=16)
ebar(xplt+0.30,np.mean(melt_mw_yr6[k_GET,120:361]),np.percentile(melt_mw_yr6[k_GET,120:361],2.5),np.percentile(melt_mw_yr6[k_GET,120:361],97.5))
plt.fill(xrect(xplt+0.25,0.2),yrect(np.mean(melt_mw_yr5[k_GET,120:361])),col5,zorder=17)
ebar(xplt+0.25,np.mean(melt_mw_yr5[k_GET,120:361]),np.percentile(melt_mw_yr5[k_GET,120:361],2.5),np.percentile(melt_mw_yr5[k_GET,120:361],97.5))

# show no BDY perturbations:
#plt.fill(xtri(xplt+0.125,0.08),ytri(np.mean(melt_mw_yr7[k_GET,120:361]),0.5),col4)
#plt.plot([xplt-0.076,xplt+0.126],[np.mean(melt_mw_yr7[k_GET,120:361]),np.mean(melt_mw_yr7[k_GET,120:361])],color='white',linewidth=2.5,zorder=4)

print([ np.mean(melt_mw_yr2[k_GET,120:361])/np.mean(melt_mw_yr1[k_GET,120:361]) , np.mean(melt_mw_yr4[k_GET,120:361])/np.mean(melt_mw_yr3[k_GET,120:361]), np.mean(melt_mw_yr6[k_GET,120:361])/np.mean(melt_mw_yr5[k_GET,120:361]) ])

print( np.mean(melt_mw_yr7[k_GET,120:361]-melt_mw_yr3[k_GET,120:361]) / (np.mean(melt_mw_yr4[k_GET,120:361]-melt_mw_yr3[k_GET,120:361]) ) )

#----------
# Dotson:
xplt=xplt+dplt
xti=np.append(xti,np.array([xplt]))

plt.fill(xrect(xplt-0.25,0.2),yrect(np.mean(melt_mw_yr2[k_DOT,120:361])),col2,zorder=21)
ebar(xplt-0.25,np.mean(melt_mw_yr2[k_DOT,120:361]),np.percentile(melt_mw_yr2[k_DOT,120:361],2.5),np.percentile(melt_mw_yr2[k_DOT,120:361],97.5))
plt.fill(xrect(xplt-0.30,0.2),yrect(np.mean(melt_mw_yr1[k_DOT,120:361])),col1,zorder=22)
ebar(xplt-0.30,np.mean(melt_mw_yr1[k_DOT,120:361]),np.percentile(melt_mw_yr1[k_DOT,120:361],2.5),np.percentile(melt_mw_yr1[k_DOT,120:361],97.5))

plt.fill(xrect(xplt+0.025,0.2),yrect(np.mean(melt_mw_yr4[k_DOT,120:361])),col4,zorder=23)
plt.plot([xplt-0.076,xplt+0.126],[np.mean(melt_mw_yr7[k_DOT,120:361]),np.mean(melt_mw_yr7[k_DOT,120:361])],color='white',linewidth=2.5,zorder=24)
ebar(xplt+0.025,np.mean(melt_mw_yr4[k_DOT,120:361]),np.percentile(melt_mw_yr4[k_DOT,120:361],2.5),np.percentile(melt_mw_yr4[k_DOT,120:361],97.5))
plt.fill(xrect(xplt-0.025,0.2),yrect(np.mean(melt_mw_yr3[k_DOT,120:361])),col3,zorder=25)
ebar(xplt-0.025,np.mean(melt_mw_yr3[k_DOT,120:361]),np.percentile(melt_mw_yr3[k_DOT,120:361],2.5),np.percentile(melt_mw_yr3[k_DOT,120:361],97.5))

plt.fill(xrect(xplt+0.30,0.2),yrect(np.mean(melt_mw_yr6[k_DOT,120:361])),col6,zorder=26)
ebar(xplt+0.30,np.mean(melt_mw_yr6[k_DOT,120:361]),np.percentile(melt_mw_yr6[k_DOT,120:361],2.5),np.percentile(melt_mw_yr6[k_DOT,120:361],97.5))
plt.fill(xrect(xplt+0.25,0.2),yrect(np.mean(melt_mw_yr5[k_DOT,120:361])),col5,zorder=27)
ebar(xplt+0.25,np.mean(melt_mw_yr5[k_DOT,120:361]),np.percentile(melt_mw_yr5[k_DOT,120:361],2.5),np.percentile(melt_mw_yr5[k_DOT,120:361],97.5))

# show no BDY perturbations:
#plt.fill(xtri(xplt+0.125,0.08),ytri(np.mean(melt_mw_yr7[k_DOT,120:361]),0.5),col4)
#plt.plot([xplt-0.076,xplt+0.126],[np.mean(melt_mw_yr7[k_DOT,120:361]),np.mean(melt_mw_yr7[k_DOT,120:361])],color='white',linewidth=2.5,zorder=24)

print([ np.mean(melt_mw_yr2[k_DOT,120:361])/np.mean(melt_mw_yr1[k_DOT,120:361]) , np.mean(melt_mw_yr4[k_DOT,120:361])/np.mean(melt_mw_yr3[k_DOT,120:361]), np.mean(melt_mw_yr6[k_DOT,120:361])/np.mean(melt_mw_yr5[k_DOT,120:361]) ])

print( np.mean(melt_mw_yr7[k_DOT,120:361]-melt_mw_yr3[k_DOT,120:361]) / (np.mean(melt_mw_yr4[k_DOT,120:361]-melt_mw_yr3[k_DOT,120:361]) ) )

#----------
# Crosson:
xplt=xplt+dplt
xti=np.append(xti,np.array([xplt]))

plt.fill(xrect(xplt-0.25,0.2),yrect(np.mean(melt_mw_yr2[k_CRO,120:361])),col2,zorder=31)
ebar(xplt-0.25,np.mean(melt_mw_yr2[k_CRO,120:361]),np.percentile(melt_mw_yr2[k_CRO,120:361],2.5),np.percentile(melt_mw_yr2[k_CRO,120:361],97.5))
plt.fill(xrect(xplt-0.30,0.2),yrect(np.mean(melt_mw_yr1[k_CRO,120:361])),col1,zorder=32)
ebar(xplt-0.30,np.mean(melt_mw_yr1[k_CRO,120:361]),np.percentile(melt_mw_yr1[k_CRO,120:361],2.5),np.percentile(melt_mw_yr1[k_CRO,120:361],97.5))

plt.fill(xrect(xplt+0.025,0.2),yrect(np.mean(melt_mw_yr4[k_CRO,120:361])),col4,zorder=33)
plt.plot([xplt-0.076,xplt+0.126],[np.mean(melt_mw_yr7[k_CRO,120:361]),np.mean(melt_mw_yr7[k_CRO,120:361])],color='white',linewidth=2.5,zorder=34)
ebar(xplt+0.025,np.mean(melt_mw_yr4[k_CRO,120:361]),np.percentile(melt_mw_yr4[k_CRO,120:361],2.5),np.percentile(melt_mw_yr4[k_CRO,120:361],97.5))
plt.fill(xrect(xplt-0.025,0.2),yrect(np.mean(melt_mw_yr3[k_CRO,120:361])),col3,zorder=35)
ebar(xplt-0.025,np.mean(melt_mw_yr3[k_CRO,120:361]),np.percentile(melt_mw_yr3[k_CRO,120:361],2.5),np.percentile(melt_mw_yr3[k_CRO,120:361],97.5))

plt.fill(xrect(xplt+0.30,0.2),yrect(np.mean(melt_mw_yr6[k_CRO,120:361])),col6,zorder=36)
ebar(xplt+0.30,np.mean(melt_mw_yr6[k_CRO,120:361]),np.percentile(melt_mw_yr6[k_CRO,120:361],2.5),np.percentile(melt_mw_yr6[k_CRO,120:361],97.5))
plt.fill(xrect(xplt+0.25,0.2),yrect(np.mean(melt_mw_yr5[k_CRO,120:361])),col5,zorder=37)
ebar(xplt+0.25,np.mean(melt_mw_yr5[k_CRO,120:361]),np.percentile(melt_mw_yr5[k_CRO,120:361],2.5),np.percentile(melt_mw_yr5[k_CRO,120:361],97.5))

# show no BDY perturbations:
#plt.fill(xtri(xplt+0.125,0.08),ytri(np.mean(melt_mw_yr7[k_CRO,120:361]),0.5),col4)
#plt.plot([xplt-0.076,xplt+0.126],[np.mean(melt_mw_yr7[k_CRO,120:361]),np.mean(melt_mw_yr7[k_CRO,120:361])],color='white',linewidth=2.5,zorder=34)

print([ np.mean(melt_mw_yr2[k_CRO,120:361])/np.mean(melt_mw_yr1[k_CRO,120:361]) , np.mean(melt_mw_yr4[k_CRO,120:361])/np.mean(melt_mw_yr3[k_CRO,120:361]), np.mean(melt_mw_yr6[k_CRO,120:361])/np.mean(melt_mw_yr5[k_CRO,120:361]) ])

print( np.mean(melt_mw_yr7[k_CRO,120:361]-melt_mw_yr3[k_CRO,120:361]) / (np.mean(melt_mw_yr4[k_CRO,120:361]-melt_mw_yr3[k_CRO,120:361]) ) )
#----------
# Thwaites:
xplt=xplt+dplt
xti=np.append(xti,np.array([xplt]))

plt.fill(xrect(xplt-0.25,0.2),yrect(np.mean(melt_mw_yr2[k_THW,120:361])),col2,zorder=41)
ebar(xplt-0.25,np.mean(melt_mw_yr2[k_THW,120:361]),np.percentile(melt_mw_yr2[k_THW,120:361],2.5),np.percentile(melt_mw_yr2[k_THW,120:361],97.5))
plt.fill(xrect(xplt-0.30,0.2),yrect(np.mean(melt_mw_yr1[k_THW,120:361])),col1,zorder=42)
ebar(xplt-0.30,np.mean(melt_mw_yr1[k_THW,120:361]),np.percentile(melt_mw_yr1[k_THW,120:361],2.5),np.percentile(melt_mw_yr1[k_THW,120:361],97.5))

plt.fill(xrect(xplt+0.025,0.2),yrect(np.mean(melt_mw_yr4[k_THW,120:361])),col4,zorder=43)
plt.plot([xplt-0.076,xplt+0.126],[np.mean(melt_mw_yr7[k_THW,120:361]),np.mean(melt_mw_yr7[k_THW,120:361])],color='white',linewidth=2.5,zorder=44)
ebar(xplt+0.025,np.mean(melt_mw_yr4[k_THW,120:361]),np.percentile(melt_mw_yr4[k_THW,120:361],2.5),np.percentile(melt_mw_yr4[k_THW,120:361],97.5))
plt.fill(xrect(xplt-0.025,0.2),yrect(np.mean(melt_mw_yr3[k_THW,120:361])),col3,zorder=45)
ebar(xplt-0.025,np.mean(melt_mw_yr3[k_THW,120:361]),np.percentile(melt_mw_yr3[k_THW,120:361],2.5),np.percentile(melt_mw_yr3[k_THW,120:361],97.5))

plt.fill(xrect(xplt+0.30,0.2),yrect(np.mean(melt_mw_yr6[k_THW,120:361])),col6,zorder=46)
ebar(xplt+0.30,np.mean(melt_mw_yr6[k_THW,120:361]),np.percentile(melt_mw_yr6[k_THW,120:361],2.5),np.percentile(melt_mw_yr6[k_THW,120:361],97.5))
plt.fill(xrect(xplt+0.25,0.2),yrect(np.mean(melt_mw_yr5[k_THW,120:361])),col5,zorder=47)
ebar(xplt+0.25,np.mean(melt_mw_yr5[k_THW,120:361]),np.percentile(melt_mw_yr5[k_THW,120:361],2.5),np.percentile(melt_mw_yr5[k_THW,120:361],97.5))

# show no BDY perturbations:
#plt.fill(xtri(xplt+0.125,0.08),ytri(np.mean(melt_mw_yr7[k_THW,120:361]),0.5),col4)
#plt.plot([xplt-0.076,xplt+0.126],[np.mean(melt_mw_yr7[k_THW,120:361]),np.mean(melt_mw_yr7[k_THW,120:361])],color='white',linewidth=2.5,zorder=44)

print([ np.mean(melt_mw_yr2[k_THW,120:361])/np.mean(melt_mw_yr1[k_THW,120:361]) , np.mean(melt_mw_yr4[k_THW,120:361])/np.mean(melt_mw_yr3[k_THW,120:361]), np.mean(melt_mw_yr6[k_THW,120:361])/np.mean(melt_mw_yr5[k_THW,120:361]) ])

print( np.mean(melt_mw_yr7[k_THW,120:361]-melt_mw_yr3[k_THW,120:361]) / (np.mean(melt_mw_yr4[k_THW,120:361]-melt_mw_yr3[k_THW,120:361]) ) )
#----------
# Pine Island:
xplt=xplt+dplt
xti=np.append(xti,np.array([xplt]))

plt.fill(xrect(xplt-0.25,0.2),yrect(np.mean(melt_mw_yr2[k_PIG,120:361])),col2,zorder=51)
ebar(xplt-0.25,np.mean(melt_mw_yr2[k_PIG,120:361]),np.percentile(melt_mw_yr2[k_PIG,120:361],2.5),np.percentile(melt_mw_yr2[k_PIG,120:361],97.5))
plt.fill(xrect(xplt-0.30,0.2),yrect(np.mean(melt_mw_yr1[k_PIG,120:361])),col1,zorder=52)
ebar(xplt-0.30,np.mean(melt_mw_yr1[k_PIG,120:361]),np.percentile(melt_mw_yr1[k_PIG,120:361],2.5),np.percentile(melt_mw_yr1[k_PIG,120:361],97.5))

plt.fill(xrect(xplt+0.025,0.2),yrect(np.mean(melt_mw_yr4[k_PIG,120:361])),col4,zorder=53)
plt.plot([xplt-0.076,xplt+0.126],[np.mean(melt_mw_yr7[k_PIG,120:361]),np.mean(melt_mw_yr7[k_PIG,120:361])],color='white',linewidth=2.5,zorder=54)
ebar(xplt+0.025,np.mean(melt_mw_yr4[k_PIG,120:361]),np.percentile(melt_mw_yr4[k_PIG,120:361],2.5),np.percentile(melt_mw_yr4[k_PIG,120:361],97.5))
plt.fill(xrect(xplt-0.025,0.2),yrect(np.mean(melt_mw_yr3[k_PIG,120:361])),col3,zorder=55)
ebar(xplt-0.025,np.mean(melt_mw_yr3[k_PIG,120:361]),np.percentile(melt_mw_yr3[k_PIG,120:361],2.5),np.percentile(melt_mw_yr3[k_PIG,120:361],97.5))

plt.fill(xrect(xplt+0.30,0.2),yrect(np.mean(melt_mw_yr6[k_PIG,120:361])),col6,zorder=56)
ebar(xplt+0.30,np.mean(melt_mw_yr6[k_PIG,120:361]),np.percentile(melt_mw_yr6[k_PIG,120:361],2.5),np.percentile(melt_mw_yr6[k_PIG,120:361],97.5))
plt.fill(xrect(xplt+0.25,0.2),yrect(np.mean(melt_mw_yr5[k_PIG,120:361])),col5,zorder=57)
ebar(xplt+0.25,np.mean(melt_mw_yr5[k_PIG,120:361]),np.percentile(melt_mw_yr5[k_PIG,120:361],2.5),np.percentile(melt_mw_yr5[k_PIG,120:361],97.5))

# show no BDY perturbations:
#plt.fill(xtri(xplt+0.125,0.08),ytri(np.mean(melt_mw_yr7[k_PIG,120:361]),0.5),col4)
#plt.plot([xplt-0.076,xplt+0.126],[np.mean(melt_mw_yr7[k_PIG,120:361]),np.mean(melt_mw_yr7[k_PIG,120:361])],color='white',linewidth=2.5,zorder=54)

print([ np.mean(melt_mw_yr2[k_PIG,120:361])/np.mean(melt_mw_yr1[k_PIG,120:361]) , np.mean(melt_mw_yr4[k_PIG,120:361])/np.mean(melt_mw_yr3[k_PIG,120:361]), np.mean(melt_mw_yr6[k_PIG,120:361])/np.mean(melt_mw_yr5[k_PIG,120:361]) ])

print( np.mean(melt_mw_yr7[k_PIG,120:361]-melt_mw_yr3[k_PIG,120:361]) / (np.mean(melt_mw_yr4[k_PIG,120:361]-melt_mw_yr3[k_PIG,120:361]) ) )
#----------
# Cosgrove:
xplt=xplt+dplt
xti=np.append(xti,np.array([xplt]))

plt.fill(xrect(xplt-0.25,0.2),yrect(np.mean(melt_mw_yr2[k_COS,120:361])),col2,zorder=61)
ebar(xplt-0.25,np.mean(melt_mw_yr2[k_COS,120:361]),np.percentile(melt_mw_yr2[k_COS,120:361],2.5),np.percentile(melt_mw_yr2[k_COS,120:361],97.5))
plt.fill(xrect(xplt-0.30,0.2),yrect(np.mean(melt_mw_yr1[k_COS,120:361])),col1,zorder=62)
ebar(xplt-0.30,np.mean(melt_mw_yr1[k_COS,120:361]),np.percentile(melt_mw_yr1[k_COS,120:361],2.5),np.percentile(melt_mw_yr1[k_COS,120:361],97.5))

plt.fill(xrect(xplt+0.025,0.2),yrect(np.mean(melt_mw_yr4[k_COS,120:361])),col4,zorder=63)
plt.plot([xplt-0.076,xplt+0.126],[np.mean(melt_mw_yr7[k_COS,120:361]),np.mean(melt_mw_yr7[k_COS,120:361])],color='white',linewidth=2.5,zorder=64)
ebar(xplt+0.025,np.mean(melt_mw_yr4[k_COS,120:361]),np.percentile(melt_mw_yr4[k_COS,120:361],2.5),np.percentile(melt_mw_yr4[k_COS,120:361],97.5))
plt.fill(xrect(xplt-0.025,0.2),yrect(np.mean(melt_mw_yr3[k_COS,120:361])),col3,zorder=65)
ebar(xplt-0.025,np.mean(melt_mw_yr3[k_COS,120:361]),np.percentile(melt_mw_yr3[k_COS,120:361],2.5),np.percentile(melt_mw_yr3[k_COS,120:361],97.5))

plt.fill(xrect(xplt+0.30,0.2),yrect(np.mean(melt_mw_yr6[k_COS,120:361])),col6,zorder=66)
ebar(xplt+0.30,np.mean(melt_mw_yr6[k_COS,120:361]),np.percentile(melt_mw_yr6[k_COS,120:361],2.5),np.percentile(melt_mw_yr6[k_COS,120:361],97.5))
plt.fill(xrect(xplt+0.25,0.2),yrect(np.mean(melt_mw_yr5[k_COS,120:361])),col5,zorder=67)
ebar(xplt+0.25,np.mean(melt_mw_yr5[k_COS,120:361]),np.percentile(melt_mw_yr5[k_COS,120:361],2.5),np.percentile(melt_mw_yr5[k_COS,120:361],97.5))

# show no BDY perturbations:
#plt.fill(xtri(xplt+0.125,0.08),ytri(np.mean(melt_mw_yr7[k_COS,120:361]),0.5),col4)
#plt.plot([xplt-0.076,xplt+0.126],[np.mean(melt_mw_yr7[k_COS,120:361]),np.mean(melt_mw_yr7[k_COS,120:361])],color='white',linewidth=2.5,zorder=64)

print([ np.mean(melt_mw_yr2[k_COS,120:361])/np.mean(melt_mw_yr1[k_COS,120:361]) , np.mean(melt_mw_yr4[k_COS,120:361])/np.mean(melt_mw_yr3[k_COS,120:361]), np.mean(melt_mw_yr6[k_COS,120:361])/np.mean(melt_mw_yr5[k_COS,120:361]) ])

print( np.mean(melt_mw_yr7[k_COS,120:361]-melt_mw_yr3[k_COS,120:361]) / (np.mean(melt_mw_yr4[k_COS,120:361]-melt_mw_yr3[k_COS,120:361]) ) )
#----------
# Abbot:
xplt=xplt+dplt
xti=np.append(xti,np.array([xplt]))

plt.fill(xrect(xplt-0.25,0.2),yrect(np.mean(melt_mw_yr2[k_ABB,120:361])),col2,zorder=71)
ebar(xplt-0.25,np.mean(melt_mw_yr2[k_ABB,120:361]),np.percentile(melt_mw_yr2[k_ABB,120:361],2.5),np.percentile(melt_mw_yr2[k_ABB,120:361],97.5))
plt.fill(xrect(xplt-0.30,0.2),yrect(np.mean(melt_mw_yr1[k_ABB,120:361])),col1,zorder=72)
ebar(xplt-0.30,np.mean(melt_mw_yr1[k_ABB,120:361]),np.percentile(melt_mw_yr1[k_ABB,120:361],2.5),np.percentile(melt_mw_yr1[k_ABB,120:361],97.5))

plt.fill(xrect(xplt+0.025,0.2),yrect(np.mean(melt_mw_yr4[k_ABB,120:361])),col4,zorder=73)
plt.plot([xplt-0.076,xplt+0.126],[np.mean(melt_mw_yr7[k_ABB,120:361]),np.mean(melt_mw_yr7[k_ABB,120:361])],color='white',linewidth=2.5,zorder=74)
ebar(xplt+0.025,np.mean(melt_mw_yr4[k_ABB,120:361]),np.percentile(melt_mw_yr4[k_ABB,120:361],2.5),np.percentile(melt_mw_yr4[k_ABB,120:361],97.5))
plt.fill(xrect(xplt-0.025,0.2),yrect(np.mean(melt_mw_yr3[k_ABB,120:361])),col3,zorder=75)
ebar(xplt-0.025,np.mean(melt_mw_yr3[k_ABB,120:361]),np.percentile(melt_mw_yr3[k_ABB,120:361],2.5),np.percentile(melt_mw_yr3[k_ABB,120:361],97.5))

plt.fill(xrect(xplt+0.30,0.2),yrect(np.mean(melt_mw_yr6[k_ABB,120:361])),col6,zorder=76)
ebar(xplt+0.30,np.mean(melt_mw_yr6[k_ABB,120:361]),np.percentile(melt_mw_yr6[k_ABB,120:361],2.5),np.percentile(melt_mw_yr6[k_ABB,120:361],97.5))
plt.fill(xrect(xplt+0.25,0.2),yrect(np.mean(melt_mw_yr5[k_ABB,120:361])),col5,zorder=77)
ebar(xplt+0.25,np.mean(melt_mw_yr5[k_ABB,120:361]),np.percentile(melt_mw_yr5[k_ABB,120:361],2.5),np.percentile(melt_mw_yr5[k_ABB,120:361],97.5))

# show no BDY perturbations:
#plt.fill(xtri(xplt+0.125,0.08),ytri(np.mean(melt_mw_yr7[k_ABB,120:361]),0.5),col4)
#plt.plot([xplt-0.076,xplt+0.126],[np.mean(melt_mw_yr7[k_ABB,120:361]),np.mean(melt_mw_yr7[k_ABB,120:361])],color='white',linewidth=2.5,zorder=74)

print([ np.mean(melt_mw_yr2[k_ABB,120:361])/np.mean(melt_mw_yr1[k_ABB,120:361]) , np.mean(melt_mw_yr4[k_ABB,120:361])/np.mean(melt_mw_yr3[k_ABB,120:361]), np.mean(melt_mw_yr6[k_ABB,120:361])/np.mean(melt_mw_yr5[k_ABB,120:361]) ])

print( np.mean(melt_mw_yr7[k_ABB,120:361]-melt_mw_yr3[k_ABB,120:361]) / (np.mean(melt_mw_yr4[k_ABB,120:361]-melt_mw_yr3[k_ABB,120:361]) ) )
#----------
# Venable:
xplt=xplt+dplt
xti=np.append(xti,np.array([xplt]))

plt.fill(xrect(xplt-0.25,0.2),yrect(np.mean(melt_mw_yr2[k_VEN,120:361])),col2,zorder=81)
ebar(xplt-0.25,np.mean(melt_mw_yr2[k_VEN,120:361]),np.percentile(melt_mw_yr2[k_VEN,120:361],2.5),np.percentile(melt_mw_yr2[k_VEN,120:361],97.5))
plt.fill(xrect(xplt-0.30,0.2),yrect(np.mean(melt_mw_yr1[k_VEN,120:361])),col1,zorder=82)
ebar(xplt-0.30,np.mean(melt_mw_yr1[k_VEN,120:361]),np.percentile(melt_mw_yr1[k_VEN,120:361],2.5),np.percentile(melt_mw_yr1[k_VEN,120:361],97.5))

plt.fill(xrect(xplt+0.025,0.2),yrect(np.mean(melt_mw_yr4[k_VEN,120:361])),col4,zorder=83)
plt.plot([xplt-0.076,xplt+0.126],[np.mean(melt_mw_yr7[k_VEN,120:361]),np.mean(melt_mw_yr7[k_VEN,120:361])],color='white',linewidth=2.5,zorder=84)
ebar(xplt+0.025,np.mean(melt_mw_yr4[k_VEN,120:361]),np.percentile(melt_mw_yr4[k_VEN,120:361],2.5),np.percentile(melt_mw_yr4[k_VEN,120:361],97.5))
plt.fill(xrect(xplt-0.025,0.2),yrect(np.mean(melt_mw_yr3[k_VEN,120:361])),col3,zorder=85)
ebar(xplt-0.025,np.mean(melt_mw_yr3[k_VEN,120:361]),np.percentile(melt_mw_yr3[k_VEN,120:361],2.5),np.percentile(melt_mw_yr3[k_VEN,120:361],97.5))

plt.fill(xrect(xplt+0.30,0.2),yrect(np.mean(melt_mw_yr6[k_VEN,120:361])),col6,zorder=86)
ebar(xplt+0.30,np.mean(melt_mw_yr6[k_VEN,120:361]),np.percentile(melt_mw_yr6[k_VEN,120:361],2.5),np.percentile(melt_mw_yr6[k_VEN,120:361],97.5))
plt.fill(xrect(xplt+0.25,0.2),yrect(np.mean(melt_mw_yr5[k_VEN,120:361])),col5,zorder=87)
ebar(xplt+0.25,np.mean(melt_mw_yr5[k_VEN,120:361]),np.percentile(melt_mw_yr5[k_VEN,120:361],2.5),np.percentile(melt_mw_yr5[k_VEN,120:361],97.5))

# show no BDY perturbations:
#plt.fill(xtri(xplt+0.125,0.08),ytri(np.mean(melt_mw_yr7[k_VEN,120:361]),0.5),col4)
#plt.plot([xplt-0.076,xplt+0.126],[np.mean(melt_mw_yr7[k_VEN,120:361]),np.mean(melt_mw_yr7[k_VEN,120:361])],color='white',linewidth=2.5,zorder=84)

print([ np.mean(melt_mw_yr2[k_VEN,120:361])/np.mean(melt_mw_yr1[k_VEN,120:361]) , np.mean(melt_mw_yr4[k_VEN,120:361])/np.mean(melt_mw_yr3[k_VEN,120:361]), np.mean(melt_mw_yr6[k_VEN,120:361])/np.mean(melt_mw_yr5[k_VEN,120:361]) ])

print( np.mean(melt_mw_yr7[k_VEN,120:361]-melt_mw_yr3[k_VEN,120:361]) / (np.mean(melt_mw_yr4[k_VEN,120:361]-melt_mw_yr3[k_VEN,120:361]) ) )
#----------
# Legend:
sq=[0.,0.,1.,1.,0.]

plt.fill(xrect(0.7*xplt,0.3),34.*np.ones((5))+sq,col1)
plt.text(0.7*xplt+0.35,34.5,'1989-2009 (A)',fontsize=8,ha='left',va='center')
plt.fill(xrect(0.7*xplt,0.3),32.*np.ones((5))+sq,col3)
plt.text(0.7*xplt+0.35,32.5,'1989-2009 (B)',fontsize=8,ha='left',va='center')
plt.fill(xrect(0.7*xplt,0.3),30.*np.ones((5))+sq,col5)
plt.text(0.7*xplt+0.35,30.5,'1989-2009 (C)',fontsize=8,ha='left',va='center')

plt.fill(xrect(0.7*xplt,0.3),27.*np.ones((5))+sq,col2)
plt.text(0.7*xplt+0.35,27.5,'2080-2100 (A)',fontsize=8,ha='left',va='center')
plt.fill(xrect(0.7*xplt,0.3),25.*np.ones((5))+sq,col4)
plt.text(0.7*xplt+0.35,25.5,'2080-2100 (B)',fontsize=8,ha='left',va='center')
plt.fill(xrect(0.7*xplt,0.3),23.*np.ones((5))+sq,col6)
plt.text(0.7*xplt+0.35,23.5,'2080-2100 (C)',fontsize=8,ha='left',va='center')

plt.text(0.5,35,'(a)',weight='bold',fontsize=8)

plt.ylim(bottom=0.0)
plt.ylabel('ice-shelf melt rate (m.w.e./yr)',size=8)
plt.yticks(fontsize=8)
ax.set_xticks(xti)
ax.set_xticklabels(('Getz','Dotson','Crosson','Thwaites','Pine Isl.','Cosgrove','Abbot','Venable'),size=8)

#plt.show()

fig.savefig('diagram_melt_ALL.jpg',dpi=300)
fig.savefig('diagram_melt_ALL.pdf')
