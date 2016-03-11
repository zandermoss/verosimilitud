#! /usr/bin/python
import Verosim as VS
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd

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
cosz=np.asarray( V.GetCosZenithEdges())



eprox_marg=np.sum(expectation,axis=1)

cosz_marg=np.sum(expectation,axis=0)
print np.sum(cosz_marg)
fig=plt.figure()

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=14)

fig,ax=plt.subplots()
"""
wid = cosz[1:] - cosz[:-1]
ax.plot(cosz[1:], cosz_marg,ls="steps",label="Expectation") 
"""
ax.set_xscale('log')
#ax.set_ylim([0,8000])
ax.set_xlabel(r'$cos(\theta_{z})$')
ax.set_ylabel(r'$E_{proxy}$')
ax.set_title("Expected Event Distribution")
#n, bins, patches = plt.hist(cosz_marg, bins=cosz, facecolor='green', alpha=0.75)
"""
#H, eprox, cosz = np.histogram2d(y, x, bins=(eprox, cosz))
#plt.set_title('pcolormesh: exact bin edges')
X, Y = np.meshgrid( cosz,eprox)
plt.pcolormesh(X, Y, expectation,norm=LogNorm(),cmap='jet')
plt.colorbar()
#ax.set_aspect('equal')
#plt.show()
"""
data=pd.read_csv("observed_events.dat",sep=" ")

data_eproxy = data['energy_proxy'].as_matrix()
data_cosz = data['cos(reconstructed_zenith)'].as_matrix()

#plt.hist(data_cosz, bins=cosz, color='Red',histtype="step",label="Data")

# Now add the legend with some customizations.
#legend = ax.legend(loc='upper left', shadow=False)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
#frame = legend.get_frame()#frame.set_facecolor('1.0')
#frame.set_facecolor('1.0')
"""
# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize('large')

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width

"""

plt.hist2d(data_eproxy, data_cosz, bins=[50,50],norm=LogNorm(),cmap="jet")
#plt.hist2d(data_eproxy, data_cosz, cmap="jet")
plt.colorbar()
plt.show()


