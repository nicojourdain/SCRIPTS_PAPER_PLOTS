import numpy as np
from scipy.stats import skewnorm
import matplotlib.pyplot as plt

mindist=999999.9

#for a in np.arange(4.94,5.10,0.02):
for a in [5.08]:

  mindist_a = 99999.9

  #for xloc in np.arange(1.96,2.12,0.02):
  for xloc in [2.02]:

    #for xscale in np.arange(1.40,1.62,0.02):
    for xscale in [1.52]:

       ECSmin = skewnorm.ppf(0.05, a, loc=xloc, scale=xscale)
       ECSmed = skewnorm.ppf(0.50, a, loc=xloc, scale=xscale)
       ECSmax = skewnorm.ppf(0.95, a, loc=xloc, scale=xscale)

       dist = np.sqrt( (ECSmin-2.0)**2 + (ECSmed-3.0)**2 + (ECSmax-5.0)**2 )

       if ( dist < mindist ):

         mindist=dist
         a_best = a
         xloc_best = xloc
         xscale_best = xscale

print(mindist, a_best, xloc_best, xscale_best)

ECS = np.arange(0,7.02,0.02)
plt.plot(ECS, skewnorm.pdf(ECS, a_best, loc=xloc_best, scale=xscale_best),color='k',linewidth=2)

ECS05=skewnorm.ppf(0.05, a_best, loc=xloc_best, scale=xscale_best); y05=skewnorm.pdf(ECS05, a_best, loc=xloc_best, scale=xscale_best)
ECS50=skewnorm.ppf(0.50, a_best, loc=xloc_best, scale=xscale_best); y50=skewnorm.pdf(ECS50, a_best, loc=xloc_best, scale=xscale_best)
ECS95=skewnorm.ppf(0.95, a_best, loc=xloc_best, scale=xscale_best); y95=skewnorm.pdf(ECS95, a_best, loc=xloc_best, scale=xscale_best)

plt.plot([ECS05,ECS05],[0,y05],color='k',linestyle='dashed')
plt.plot([ECS50,ECS50],[0,y50],color='k',linestyle='dashed')
plt.plot([ECS95,ECS95],[0,y95],color='k',linestyle='dashed')

plt.xlabel('ECS (°C)')
plt.ylabel('Probability')
plt.xlim([0,7])
plt.ylim([0,0.5])

# values for our CMIP6 models:
GECS=np.array([ 4.7, 3.9, 5.6, 5.2, 4.8, 4.8, 4.8, 3.9, 2.6, 3.1, 1.9, 4.6, 3.0, 3.2, 2.5, 5.3 ])
PDF=np.zeros(np.shape(GECS))
for kk in np.arange(np.size(GECS)):
  PDF[kk]=skewnorm.pdf(GECS[kk], a_best, loc=xloc_best, scale=xscale_best)
  print(GECS[kk],PDF[kk])
print('percentiles ECS = ',np.percentile(GECS,5), np.percentile(GECS,50), np.percentile(GECS,95))
print('mean ECS = ',np.mean(GECS))
print('mean weighted ECS = ',np.sum(GECS*PDF)/np.sum(PDF))

# estimate percentiles of the weighted distribution:
distrib=np.array([])
for kk in np.arange(np.size(GECS)):
  for kd in np.arange(np.round(PDF[kk]*100)):
    distrib=np.append(distrib,GECS[kk])
print('percentiles weighted ECS = ',np.percentile(distrib,5), np.percentile(distrib,50), np.percentile(distrib,95))
print('mean weighted ECS (other method) = ',np.mean(distrib))

def triangle(val):
   w=0.1
   x=[val-w,val,val+w,val-w]
   y=[0.5,0.475,0.5,0.5]
   return[x,y]

plt.fill(triangle(np.percentile(GECS,5))[0],triangle(np.percentile(GECS,5))[1],color='orange')
plt.fill(triangle(np.percentile(GECS,50))[0],triangle(np.percentile(GECS,50))[1],color='orange')
plt.fill(triangle(np.percentile(GECS,95))[0],triangle(np.percentile(GECS,95))[1],color='orange')

plt.fill(triangle(np.percentile(distrib,5))[0],triangle(np.percentile(distrib,5))[1],color='blue')
plt.fill(triangle(np.percentile(distrib,50))[0],triangle(np.percentile(distrib,50))[1],color='blue')
plt.fill(triangle(np.percentile(distrib,95))[0],triangle(np.percentile(distrib,95))[1],color='blue')

plt.savefig("skew_normal_distrib.pdf", format="pdf")

#plt.show()
