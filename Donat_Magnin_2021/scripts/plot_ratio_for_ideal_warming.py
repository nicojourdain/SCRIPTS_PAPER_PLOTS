import numpy as np
import matplotlib
from matplotlib import pyplot as plt

lab=[ 'ABB', 'COS', 'PIG', 'THW', 'CRO', 'DOT', 'GET' ]
lab2=[ 'Abbot', 'Corgrove', 'PIG', 'Thwaites', 'Crosson', 'Dotson', 'Getz' ]
col=[ 'maroon', 'tab:red', 'tab:orange', 'gold', 'aqua', 'tab:cyan', 'darkblue' ]
base_snf=np.array([943., 372., 521., 989., 1339., 830., 978.])  # for 2080-2100 in mm.w.e/yr
base_mlt=np.array([577., 588., 455., 244.,  183., 292., 333.])   # for 2080-2100 in mm.w.e/yr
Nisf=np.size(base_snf)

pres_snf=np.array([790., 290., 407., 809., 1055., 669., 786.])
pres_mlt=np.array([ 54.,  80.,  79.,  26.,   17.,  21.,  22.])

thresh1=0.60
thresh2=0.70
thresh3=0.85

warming=np.arange(-3.6,6.1,0.1) # warming compared to 2080-2100

snf=np.zeros((Nisf,np.size(warming)))
mlt=np.zeros((Nisf,np.size(warming)))
ratio=np.zeros((Nisf,np.size(warming)))
dT1=np.zeros((Nisf))
dT2=np.zeros((Nisf))

# ratio from Clausius-Clapeyron and Trusel (2015):
for kisf in np.arange(Nisf):
  #snf[kisf,:] = base_snf[kisf]*(1.0+0.095*warming)
  snf[kisf,:] = base_snf[kisf]*np.exp(17.65*warming/243.04)
  mlt[kisf,:] = base_mlt[kisf]*np.exp(0.4557*warming) 
  ratio[kisf,:] = mlt[kisf,:]/snf[kisf,:]
  k1 = np.argmin( np.abs( ratio[kisf,:] - thresh1 ) )
  k2 = np.argmin( np.abs( ratio[kisf,:] - thresh3 ) )
  dT1[kisf] = warming[k1]+3.6
  dT2[kisf] = warming[k2]+3.6

# calculate the ratio tho match all present and future points as best as possible:
# ratio(kisf,deltaT)=ration(kisf,0)*exp(AA*deltaT)
AA = np.mean( np.log(base_mlt*pres_snf/(base_snf*pres_mlt)) ) / 3.6
print 'AA=', AA
ratio2=np.zeros(np.shape(ratio))
for kisf in np.arange(Nisf):
  # expression so that the projected simulation matches (as for ratio):
  ratio2[kisf,:] = base_mlt[kisf]/base_snf[kisf] * np.exp(AA*warming)
  k1 = np.argmin( np.abs( ratio2[kisf,:] - thresh1 ) )
  k2 = np.argmin( np.abs( ratio2[kisf,:] - thresh3 ) )
  dT1[kisf] = np.minimum(dT1[kisf],warming[k1]+3.6)
  dT2[kisf] = np.maximum(dT2[kisf],warming[k2]+3.6)
  #lab[kisf] = lab[kisf]+'  ('+np.round(dT1[kisf],1).astype('str')+'-'+np.round(dT2[kisf],1).astype('str')+'$^\circ$C)'

#====================
fig, ax = plt.subplots()

shift=np.array([0.05, 0.0, 0.0, 0.1, 0.0, 0.0, 0.05])
for kisf in np.arange(Nisf):
  # ratio from Clausius-Clapeyron and Trusel (2015):
  ax.plot(warming+3.6,ratio[kisf,:],'-', label=lab2[kisf],linewidth=1.1,color=col[kisf],zorder=30+kisf)  # shift warming reference to present-day
  ax.plot(warming+3.6,ratio2[kisf,:],'--',linewidth=0.8,color=col[kisf],zorder=20+kisf)  # shift warming reference to present-day
  ax.fill([dT1[kisf],dT2[kisf],dT2[kisf],dT1[kisf],dT1[kisf]],[-0.1,-0.1,0.0,0.0,-0.1]+shift[kisf],col[kisf])

ax.set_xticks(np.arange(0,10,1.0))
ax.set_yticks(np.arange(0,3.2,0.2))
plt.xlim(0,9.0)
plt.ylim(-0.1,3.1)
ax.legend(fontsize=8).set_zorder(100)
ax.plot([0,9],[thresh1,thresh1],'--',linewidth=0.6,color='black',zorder=1)
ax.plot([0,9],[thresh2,thresh2],'-',linewidth=0.6,color='black',zorder=0)
ax.plot([0,9],[thresh3,thresh3],'--',linewidth=0.6,color='black',zorder=2)

# warming scenarios (see Table 12.2 in IPCC-AR5-ch.12):
ax.plot([3.6,3.6],[-0.1,3.1],'-',linewidth=0.6,color='gray',zorder=3)
ax.text(3.6,3.11,'2090 rcp8.5',fontsize=7,color='gray',rotation=45.0,horizontalalignment='left',verticalalignment='bottom')
ax.plot([1.9,1.9],[-0.1,3.1],':',linewidth=0.6,color='gray',zorder=4)
ax.text(1.9,3.11,'2055 rcp8.5',fontsize=7,color='gray',rotation=45.0,horizontalalignment='left',verticalalignment='bottom')
ax.plot([6.4,6.4],[-0.1,3.1],':',linewidth=0.6,color='gray',zorder=5)
ax.text(6.4,3.11,'2190 rcp8.5',fontsize=7,color='gray',rotation=45.0,horizontalalignment='left',verticalalignment='bottom')
ax.plot([7.7,7.7],[-0.1,3.1],':',linewidth=0.6,color='gray',zorder=6)
ax.text(7.7,3.11,'2290 rcp8.5',fontsize=7,color='gray',rotation=45.0,horizontalalignment='left',verticalalignment='bottom')
ax.plot([0.9,0.9],[-0.1,3.1],':',linewidth=0.6,color='gray',zorder=7)
ax.text(0.9,3.11,'2090 rcp2.6',fontsize=7,color='gray',rotation=45.0,horizontalalignment='left',verticalalignment='bottom')
ax.plot([1.7,1.7],[-0.1,3.1],':',linewidth=0.6,color='gray',zorder=8)
ax.text(1.7,3.11,'2090 rcp4.5',fontsize=7,color='gray',rotation=45.0,horizontalalignment='left',verticalalignment='bottom')
# CMIP6 SSP scenarios
ax.plot([5.056,5.056],[-0.1,3.1],'--',linewidth=0.4,color='rosybrown',zorder=10)
ax.text(5.056,3.11,'2090 ssp585',fontsize=7,color='rosybrown',rotation=45.0,horizontalalignment='left',verticalalignment='bottom')
ax.plot([2.825,2.825],[-0.1,3.1],'--',linewidth=0.4,color='rosybrown',zorder=11)
ax.text(2.825,3.11,'2090 ssp245',fontsize=7,color='rosybrown',rotation=45.0,horizontalalignment='left',verticalalignment='bottom')
ax.plot([1.713,1.713],[-0.1,3.1],'--',linewidth=0.4,color='rosybrown',zorder=12)
ax.text(1.713-0.13,3.15,'2090 ssp126',fontsize=7,color='rosybrown',rotation=45.0,horizontalalignment='left',verticalalignment='bottom')

ax.set_xlabel('warming from present-day ($^\circ$C)',size=12)
ax.set_ylabel('melt/snowfall',size=12)
ax.tick_params(axis='both', which='major', labelsize=10)

fig.savefig('mlt_to_snowfall_ratio_ideal_warming.pdf')
