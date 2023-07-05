import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

reg=[ 'ABB', 'COS', 'PIG', 'THW', 'CRO', 'DOT', 'GET' ]
reg2=[ 'Abbot', 'Corgrove', 'PIG', 'Thwaites', 'Crosson', 'Dotson', 'Getz' ]
linthic=[ 0.9, 0.7, 0.5, 0.9, 0.7, 0.5, 0.9 ]
linstyl=[ ':', '--', '-', ':', '--', '-', ':' ]
Nreg=np.size(reg)

Nyrspin=12

runoff=np.zeros((Nreg,Nyrspin,230,276))

fig, ax = plt.subplots()

spindur=np.arange(1,Nyrspin+1,1)

for kreg in np.arange(Nreg):

  file_msk='msk_isf_'+reg[kreg]+'.npy'
  msk_isf=np.load(file_msk)

  for kspin in spindur:
    file_MAR='DJF1998_'+kspin.astype('str')+'yr_spinup.nc'
    ncS=xr.open_dataset(file_MAR)
    rof=ncS['rof'].values[0,0,:,:]
    runoff[kreg,kspin-1,:,:]=rof*msk_isf

  rospin=np.nansum(np.nansum(runoff[kreg,:,:,:],axis=2),axis=1)/np.nansum(np.nansum(msk_isf,axis=1),axis=0) 
  #rospin=np.nansum(np.nansum(runoff[kreg,:,:,:],axis=2),axis=1)
  ax.plot(spindur,rospin,linstyl[kreg],label=reg2[kreg],linewidth=linthic[kreg])

ax.legend()
ax.set_xlabel('spin-up duration (yr)',size=10)
ax.set_ylabel('DJF Net Liquid (mm.w.eq/day)',size=10)
ax.tick_params(axis='both', which='major', labelsize=10)
fig.savefig('spin_up_DJF1998.pdf')
