import numpy as np
import matplotlib.pyplot as plt

##################################''
##################################''
##################################''

fig, axs = plt.subplots(nrows=3,ncols=2,figsize=(18.0,21.0),gridspec_kw={'height_ratios': [1, 1, 0.333333333]})
axs = axs.ravel()

modref=['CESM2', 'CNRM-CM6-1', 'MPI-ESM1-2-HR', 'UKESM1-0-LL', 'ACCESS1-3']
abbrev=['CES'  , 'CNR'       , 'MPI'          , 'UKE'        , 'ACC'      ]
Nref=np.size(modref)

Nmod=7

SMBref_grn=np.zeros((Nref))
SMBref_isf=np.zeros((Nref))
MLTref_grn=np.zeros((Nref))
MLTref_isf=np.zeros((Nref))
SMBmod_grn=np.zeros((Nref,Nmod))
SMBmod_isf=np.zeros((Nref,Nmod))
MLTmod_grn=np.zeros((Nref,Nmod))
MLTmod_isf=np.zeros((Nref,Nmod))

for kref in np.arange(Nref):

  filein='stats_other_models_to_'+modref[kref]+'.npz'
  ff = np.load(filein)

  # a- SMB over grounded ice sheet
  SMBref_grn[kref]=ff['mean_smb_grn_2015_2100']
  SMBmod_grn[kref,:]=ff['bias_smb_grn_2015_2100']+ff['mean_smb_grn_2015_2100']

  # b- SMB over ice shelves
  SMBref_isf[kref]=ff['mean_smb_isf_2015_2100']
  SMBmod_isf[kref,:]=ff['bias_smb_isf_2015_2100']+ff['mean_smb_isf_2015_2100']

  # c- MELT over grounded ice sheet
  MLTref_grn[kref]=ff['mean_me_grn_2015_2100']
  MLTmod_grn[kref,:]=ff['bias_me_grn_2015_2100']+ff['mean_me_grn_2015_2100']

  # d- MELT over ice shelves
  MLTref_isf[kref]=ff['mean_me_isf_2015_2100']
  MLTmod_isf[kref,:]=ff['bias_me_isf_2015_2100']+ff['mean_me_isf_2015_2100']

  col=ff['colorA']
  mod=ff['modelA']

  ff.close()

#--------------------------------------------------------------------
# prepare polygons for plots :

def data_to_polygon(data):
   """ Function to extract the coordinates of a polygon from a data vector"""
   Nd=np.size(data)
   theta=np.arange(0.,2*np.pi,2*np.pi/Nd)+np.pi/2.
   xd=data*np.cos(theta)
   yd=data*np.sin(theta)
   # to close polygon:
   xd=np.append(xd,xd[0])
   yd=np.append(yd,yd[0])
   return([xd,yd])


sc0 = np.round(np.max(SMBref_grn)/1.e2)*1.e2
sc1 = np.round(np.max(SMBref_isf*(-1))/10.)*10.
sc2 = np.round(np.max(MLTref_grn)/1.e2)*1.e2
sc3 = np.round(np.max(MLTref_isf)/1.e2)*1.e2
[x0,y0] = data_to_polygon( np.ones(np.shape(SMBref_grn)) * sc0 )
[x1,y1] = data_to_polygon( np.ones(np.shape(SMBref_isf)) * sc1 )
[x2,y2] = data_to_polygon( np.ones(np.shape(MLTref_grn)) * sc2 )
[x3,y3] = data_to_polygon( np.ones(np.shape(MLTref_isf)) * sc3 )

[xSMBref_grn,ySMBref_grn] = data_to_polygon(SMBref_grn)
[xSMBref_isf,ySMBref_isf] = data_to_polygon(np.maximum(SMBref_isf*(-1),np.zeros(np.shape(SMBref_isf))))
[xMLTref_grn,yMLTref_grn] = data_to_polygon(MLTref_grn)
[xMLTref_isf,yMLTref_isf] = data_to_polygon(MLTref_isf)

#--------------------------------------------------------------------
# plot :

print(SMBref_isf)

# a- SMB over grounded ice sheet
fac=1.4
for kref in np.arange(Nref):
  axs[0].plot([0.,fac*x0[kref]],[0.,fac*y0[kref]],color='grey',linewidth=0.5,dashes=[4, 4],zorder=kref+1)
axs[0].plot(x0*0.25*fac,y0*0.25*fac,color='k',linewidth=0.5,zorder=Nref+2)
axs[0].plot(x0*0.50*fac,y0*0.50*fac,color='k',linewidth=0.5,zorder=Nref+3)
axs[0].plot(x0*0.75*fac,y0*0.75*fac,color='k',linewidth=0.5,zorder=Nref+4)
axs[0].plot(x0*1.00*fac,y0*1.00*fac,color='k',linewidth=0.5,zorder=Nref+5)
#-
val=np.sqrt(x0[2]**2+y0[2]**2)*fac
axs[0].text(0.,y0[2]*fac*0.25,round(val*0.25),color='grey',fontsize=12,ha='center',va='center',backgroundcolor='white',zorder=20)
axs[0].text(0.,y0[2]*fac*0.50,round(val*0.50),color='grey',fontsize=12,ha='center',va='center',backgroundcolor='white',zorder=21)
axs[0].text(0.,y0[2]*fac*0.75,round(val*0.75),color='grey',fontsize=12,ha='center',va='center',backgroundcolor='white',zorder=22)
axs[0].text(0.,y0[2]*fac*1.00,round(val*1.00),color='grey',fontsize=12,ha='center',va='center',backgroundcolor='white',zorder=23)
#-
for kref in np.arange(Nref):
  axs[0].text(x0[kref]*fac*1.10,y0[kref]*fac*1.05,abbrev[kref],color='k',ha='center',va='center',zorder=50+kref,fontsize=16)
#-
axs[0].plot(xSMBref_grn,ySMBref_grn,color='k',linewidth=4,marker='o',markersize=20,zorder=30)
for kmod in np.arange(Nmod):
  [xx,yy] = data_to_polygon(SMBmod_grn[:,kmod])
  axs[0].plot(xx,yy,color=col[kmod],linewidth=2,marker='o',markersize=10,zorder=31+kmod)
#-
axs[0].set_aspect('equal')
axs[0].set_axis_off()
axs[0].set_title('(a) SMB anomaly over the grounded ice (Gt/yr)',fontsize=16,fontweight='bold')

# b- SMB over ice shelves
fac=2.0
for kref in np.arange(Nref):
  axs[1].plot([0.,fac*x1[kref]],[0.,fac*y1[kref]],color='grey',linewidth=0.5,dashes=[4, 4],zorder=kref+1)
axs[1].plot(x1*0.25*fac,y1*0.25*fac,color='k',linewidth=0.5,zorder=Nref+2)
axs[1].plot(x1*0.50*fac,y1*0.50*fac,color='k',linewidth=0.5,zorder=Nref+3)
axs[1].plot(x1*0.75*fac,y1*0.75*fac,color='k',linewidth=0.5,zorder=Nref+4)
axs[1].plot(x1*1.00*fac,y1*1.00*fac,color='k',linewidth=0.5,zorder=Nref+5)
#-
val=-np.sqrt(x1[2]**2+y1[2]**2)*fac
axs[1].text(0.,y1[2]*fac*0.25,round(val*0.25),color='grey',fontsize=12,ha='center',va='center',backgroundcolor='white',zorder=20)
axs[1].text(0.,y1[2]*fac*0.50,round(val*0.50),color='grey',fontsize=12,ha='center',va='center',backgroundcolor='white',zorder=21)
axs[1].text(0.,y1[2]*fac*0.75,round(val*0.75),color='grey',fontsize=12,ha='center',va='center',backgroundcolor='white',zorder=22)
axs[1].text(0.,y1[2]*fac*1.00,round(val*1.00),color='grey',fontsize=12,ha='center',va='center',backgroundcolor='white',zorder=23)
#-
for kref in np.arange(Nref):
  axs[1].text(x1[kref]*fac*1.10,y1[kref]*fac*1.05,abbrev[kref],color='k',ha='center',va='center',zorder=50+kref,fontsize=16)
#-
axs[1].plot(xSMBref_isf,ySMBref_isf,color='k',linewidth=4,marker='o',markersize=20,zorder=30)
for kmod in np.arange(Nmod):
  [xx,yy] = data_to_polygon(np.maximum(SMBmod_isf[:,kmod]*(-1),np.zeros(np.shape(SMBmod_isf[:,kmod]))))
  axs[1].plot(xx,yy,color=col[kmod],linewidth=2,marker='o',markersize=10,zorder=31+kmod)
#-
axs[1].set_aspect('equal')
axs[1].set_axis_off()
axs[1].set_title('(b) SMB anomaly over the ice shelves (Gt/yr)',fontsize=16,fontweight='bold')

# c- MELT over grounded ice sheet
fac=1.4
for kref in np.arange(Nref):
  axs[2].plot([0.,fac*x2[kref]],[0.,fac*y2[kref]],color='grey',linewidth=0.5,dashes=[4, 4],zorder=kref+1)
axs[2].plot(x2*0.25*fac,y2*0.25*fac,color='k',linewidth=0.5,zorder=Nref+2)
axs[2].plot(x2*0.50*fac,y2*0.50*fac,color='k',linewidth=0.5,zorder=Nref+3)
axs[2].plot(x2*0.75*fac,y2*0.75*fac,color='k',linewidth=0.5,zorder=Nref+4)
axs[2].plot(x2*1.00*fac,y2*1.00*fac,color='k',linewidth=0.5,zorder=Nref+5)
#-
val=np.sqrt(x2[2]**2+y2[2]**2)*fac
axs[2].text(0.,y2[2]*fac*0.25,round(val*0.25),color='grey',fontsize=12,ha='center',va='center',backgroundcolor='white',zorder=20)
axs[2].text(0.,y2[2]*fac*0.50,round(val*0.50),color='grey',fontsize=12,ha='center',va='center',backgroundcolor='white',zorder=21)
axs[2].text(0.,y2[2]*fac*0.75,round(val*0.75),color='grey',fontsize=12,ha='center',va='center',backgroundcolor='white',zorder=22)
axs[2].text(0.,y2[2]*fac*1.00,round(val*1.00),color='grey',fontsize=12,ha='center',va='center',backgroundcolor='white',zorder=23)
#-
for kref in np.arange(Nref):
  axs[2].text(x2[kref]*fac*1.10,y2[kref]*fac*1.05,abbrev[kref],color='k',ha='center',va='center',zorder=50+kref,fontsize=16)
#-
axs[2].plot(xMLTref_grn,yMLTref_grn,color='k',linewidth=4,marker='o',markersize=20,zorder=30)
for kmod in np.arange(Nmod):
  [xx,yy] = data_to_polygon(MLTmod_grn[:,kmod])
  axs[2].plot(xx,yy,color=col[kmod],linewidth=2,marker='o',markersize=10,zorder=31+kmod)
#-
axs[2].set_aspect('equal')
axs[2].set_axis_off()
axs[2].set_title('(c) Surface melting anomaly over the grounded ice (Gt/yr)',fontsize=16,fontweight='bold')

# d- MELT over ice shelves
fac=1.5
for kref in np.arange(Nref):
  axs[3].plot([0.,fac*x3[kref]],[0.,fac*y3[kref]],color='grey',linewidth=0.5,dashes=[4, 4],zorder=kref+1)
axs[3].plot(x3*0.25*fac,y3*0.25*fac,color='k',linewidth=0.5,zorder=Nref+2)
axs[3].plot(x3*0.50*fac,y3*0.50*fac,color='k',linewidth=0.5,zorder=Nref+3)
axs[3].plot(x3*0.75*fac,y3*0.75*fac,color='k',linewidth=0.5,zorder=Nref+4)
axs[3].plot(x3*1.00*fac,y3*1.00*fac,color='k',linewidth=0.5,zorder=Nref+5)
#-
val=np.sqrt(x3[2]**2+y3[2]**2)*fac
axs[3].text(0.,y3[2]*fac*0.25,round(val*0.25),color='grey',fontsize=12,ha='center',va='center',backgroundcolor='white',zorder=20)
axs[3].text(0.,y3[2]*fac*0.50,round(val*0.50),color='grey',fontsize=12,ha='center',va='center',backgroundcolor='white',zorder=21)
axs[3].text(0.,y3[2]*fac*0.75,round(val*0.75),color='grey',fontsize=12,ha='center',va='center',backgroundcolor='white',zorder=22)
axs[3].text(0.,y3[2]*fac*1.00,round(val*1.00),color='grey',fontsize=12,ha='center',va='center',backgroundcolor='white',zorder=23)
#-
for kref in np.arange(Nref):
  axs[3].text(x3[kref]*fac*1.10,y3[kref]*fac*1.05,abbrev[kref],color='k',ha='center',va='center',zorder=50+kref,fontsize=16)
#-
axs[3].plot(xMLTref_isf,yMLTref_isf,color='k',linewidth=4,marker='o',markersize=20,zorder=30)
for kmod in np.arange(Nmod):
  [xx,yy] = data_to_polygon(MLTmod_isf[:,kmod])
  axs[3].plot(xx,yy,color=col[kmod],linewidth=2,marker='o',markersize=10,zorder=31+kmod)
#-
axs[3].set_aspect('equal')
axs[3].set_axis_off()
axs[3].set_title('(d) Surface melting anomaly over the ice shelves (Gt/yr)',fontsize=16,fontweight='bold')

#- homemade legend:
axs[4].plot([0.2,0.35],[1.,1.],color='k',linewidth=4,marker='o',markersize=20)
axs[4].text(0.4,1,'Original MAR simulation driven by:',color='k',ha='left',va='center',fontsize=16,fontweight='bold')
for kref in np.arange(Nref):
  axs[4].text(0.5,1-(kref+1)*0.15,'- '+modref[kref]+' ('+abbrev[kref]+')',color='k',ha='left',va='center',fontsize=15)
axs[4].set_xlim(0,1.1)
axs[4].set_ylim(-0.1,1.1)
axs[4].set_axis_off()

axs[5].text(0,1,'Reconstructions of MAR simulation driven by:',color='k',ha='left',va='center',fontsize=16,fontweight='bold')
for kmod in np.arange(Nmod):
  axs[5].plot([0.1,0.25],[1-(kmod+1)*0.15,1-(kmod+1)*0.15],color=col[kmod],linewidth=2,marker='o',markersize=10)
  axs[5].text(0.3,1-(kmod+1)*0.15,mod[kmod],color='k',ha='left',va='center',fontsize=15)
axs[5].set_xlim(0,1.1)
axs[5].set_ylim(-0.1,1.1)
axs[5].set_axis_off()

#------------------------------------------------------------------------
figname='polygons_other_models_to_all.pdf'

fig.savefig(figname)
