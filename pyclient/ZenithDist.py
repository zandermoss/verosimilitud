#! /usr/bin/python
import Verosim as VS
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd
from numpy import sqrt

V = VS.Verosim(4)



raw_expectation= V.Likelihood()
dims=V.GetExpDims()

#deserialze

datasize=1
strides=[]
for i in range(0,2):
	strides.append(datasize)
	datasize*=dims[2-(i+1)]

index=0

expectation = np.zeros(dims)

def index(indices,raw_expectation):
	myindex=0
	for i in range(0,2):
		myindex+=strides[i]*indices[2-(i+1)]
	
	return raw_expectation[myindex] 
olv=0
for x in range(0,dims[0]):
	for y in range(0,dims[1]):
		myind=[x,y]		
		expectation[x,y]=index(myind,raw_expectation)
		olv+= expectation[x,y]
#print expectation


eprox= np.asarray(V.GetEproxEdges())
#cosz=np.asarray( V.GetCosZenithEdges())
cosz= [-1. , -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,  0.,   0.1]

print "COSZ: ", cosz


fig,ax=plt.subplots()
data=np.genfromtxt("./observed_events.dat")
data_eproxy = data[:,2]
data_cosz = data[:,3]

bins_center_costh = [ (cosz[i]+cosz[i+1])/2.0 for i in range(len(cosz)-1) ]
print data_cosz
print len(data_cosz)
print "COSZ2: ", cosz
n,b,p=ax.hist(data_cosz, bins=cosz, color='Green',histtype="step",label="Data",lw=1)

print n

print "SUM: ",np.sum(n)

error = lambda x : sqrt(x)
ax.errorbar(bins_center_costh, n, yerr=error(n), color = "black",linewidth=0,elinewidth=2)


ax.plot(bins_center_costh,n,'*')



eprox_marg=np.sum(expectation,axis=1)

cosz_marg=np.sum(expectation,axis=0)
print np.sum(cosz_marg)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=14)



ax.hist(bins_center_costh,bins = cosz, weights = cosz_marg,histtype = 'step',color = "Red", lw =2, alpha=1.0,label = 'Expectation')


weaver=np.genfromtxt("./weaver.dat")
weaver_exp = weaver[:,0]
weaver_data = weaver[:,1]


print "WLEN: ",len(weaver_exp)
print "CLEN: ",len(bins_center_costh)

ax.hist(bins_center_costh,bins = cosz, weights = weaver_exp,histtype = 'step',color = "Red", lw =2, alpha=1.0,label = 'Weaver Expectation',ls='dashed')

ax.hist(bins_center_costh,bins = cosz, weights = weaver_data,histtype = 'step',color = "Green", lw =2, alpha=1.0,label = 'Weaver Data',ls='dashed')


# Now add the legend with some customizations.
legend = ax.legend(loc='upper left', shadow=False)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame = legend.get_frame()
frame.set_facecolor('1.0')

# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize('large')

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width


#ax.set_yscale('log')
ax.set_ylim([0,8000])
ax.set_xlabel(r'$cos(\theta_{z})$')
ax.set_ylabel('Events')
ax.set_title("Event Distribution over Cos (Zenith Angle)")
#plt.semilogy()
plt.show()
#H, eprox, cosz = plt.histogram2d(data_eproxy, data_cosz, bins=(eprox, cosz))
#plt.show()


