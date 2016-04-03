#! /usr/bin/python
import Verosim as VS
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd

V = VS.Verosim(4)

cuts=np.zeros(2)
cuts[0]=6
cuts[1]=50-16

V.SetEproxCuts(cuts)

V.CalculateExpectation()

init=np.zeros(2)
init[0]=1
init[1]=0

retvec=V.Chi2MinNuisance(init)
print retvec

nuisance=np.zeros(1)
nuisance[0]=retvec[1]
#nuisance[1]=retvec[2]

eprox=np.asarray( V.GetEproxEdges())
cosz=np.asarray(V.GetCosZenithEdges())

#raw_expectation= V.GetExpectationVec()
raw_expectation= V.GetExpectationVec(nuisance)
#raw_data= V.GetDataVec()
raw_data= V.GetDataVec(1.5)
dims=V.GetExpDims()
print "DIMS: ", dims

#deserialze
datasize=1
strides=[]
for i in range(0,2):
	strides.append(datasize)
	datasize*=dims[2-(i+1)]

index=0

expectation = np.zeros(dims)
data = np.zeros(dims)

def index(indices,raw_expectation):
	myindex=0
	for i in range(0,2):
		myindex+=strides[i]*indices[2-(i+1)]
	#print raw_expectation[myindex]	
	return raw_expectation[myindex] 

olv=0

for x in range(0,dims[0]):
	for y in range(0,dims[1]):
		myind=[x,y]		
		expectation[x,y]=index(myind,raw_expectation)
		data[x,y]=index(myind,raw_data)
		olv+= expectation[x,y]

nopert=np.load("nopert.npz")
nopert_eprox_exp_marg=np.sum(nopert['exp'],axis=1)
nopert_cosz_exp_marg=np.sum(nopert['exp'],axis=0)
nopert_eprox_dat_marg=np.sum(nopert['dat'],axis=1)
nopert_cosz_dat_marg=np.sum(nopert['dat'],axis=0)

eprox_marg=np.sum(expectation,axis=1)
cosz_marg=np.sum(expectation,axis=0)
#print cosz_marg

data_eprox_marg=np.sum(data,axis=1)
data_cosz_marg=np.sum(data,axis=0)
#print data_cosz_marg


#print np.sum(cosz_marg)
#print np.sum(data_cosz_marg)
fig=plt.figure()

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=14)


#print np.divide(cosz_marg,data_cosz_marg)

fig,ax=plt.subplots()
#wid = eprox[1:] - eprox[:-1]
ax.plot(cosz[1:], np.divide(cosz_marg,data_cosz_marg),ls="steps",label="Expectation/Data Pert")
ax.plot(cosz[1:], np.divide(nopert_cosz_exp_marg,nopert_cosz_dat_marg),ls="steps",label="Expectation/Data Nopert",color="red")
#ax.plot(eprox[1:], nopert_eprox_dat_marg,ls="steps",label="Nopert")
#ax.plot(eprox[1:], data_eprox_marg,ls="steps",label="Pert")
#ax.set_xscale('log')
ax.set_ylim([0.8,1.2])
ax.set_xlabel("cosz")
ax.set_ylabel('Ratio')
ax.set_title("Ratio in cosz")
#n, bins, patches = plt.hist(cosz_marg, bins=cosz, facecolor='green', alpha=0.75)

"""
#H, eprox, cosz = np.histogram2d(y, x, bins=(eprox, cosz))
fig=plt.figure()
#plt.set_title('pcolormesh: exact bin edges')
X, Y = np.meshgrid( cosz,eprox)
plt.pcolormesh(X, Y, expectation,norm=LogNorm(),cmap='PuBu_r')
plt.colorbar()
#ax.set_aspect('equal')
"""
#plt.show()



# Now add the legend with some customizations.
legend = ax.legend(loc='upper right', shadow=False)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame = legend.get_frame()
frame.set_facecolor('1.0')

# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize('large')

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width

#plt.semilogy()

plt.show()
#H, eprox, cosz = plt.histogram2d(data_eproxy, data_cosz, bins=(eprox, cosz))
#plt.show()


