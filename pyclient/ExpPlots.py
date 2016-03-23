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

np.savez("exp_and_edge",exp=expectation,ep=eprox,cos=cosz)

#eprox_marg=np.sum(expectation,axis=1)

