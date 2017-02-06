#! /usr/bin/python

import PhysConst as pc
import simple_propagate as sp
import numpy as np
import matplotlib.pyplot as plt
import PMNSGen
import random
from math import pi
import math
from tqdm import tqdm

import Verosim

#Set paths to data, fluxes, and responses.
data_path="/home/pinkpig/physics/neutrino_decay/nuisance_runs/neutrino_decay/verosimilitud/data/IC86SterileNeutrinoDataRelease/data/observed_events.dat"
flux_path="/home/pinkpig/physics/neutrino_decay/nuisance_runs/neutrino_decay/verosimilitud/data/IC86SterileNeutrinoDataRelease/atmospheric_flux/averaged/PolyGonato_QGSJET-II-04.h5"
effective_area_path="/home/pinkpig/physics/neutrino_decay/nuisance_runs/neutrino_decay/verosimilitud/data/IC86SterileNeutrinoDataRelease/systematics_response_arrays/"


#Load oscillation probability arrays
osc_file = np.load("3flav_matter_osc.npz")
nu_prob = osc_file['nu']
nubar_prob = osc_file['nubar']
#nu_prob = np.ones(len(nu_prob))
#nubar_prob = np.ones(len(nubar_prob))

#Calling Verosimilitud
V=Verosim.Verosim(4,data_path, flux_path, effective_area_path, nu_prob, nubar_prob)

#Testing Chi2
"""
  There are 5 nuisance parameters at the moment.
  #  Name                        Mean    Sigma
  --------------------------------------------
  1) Normalization coefficient   1.0     0.4
  2) Spectral index              0.0     0.05
  3) Kaon-pion ratio             1.0     0.1
  4) Nu-nubar ratio              1.0     0.025
  5) DOM efficiency              1.0     0.3
"""

#Testing chi2 estimation at mean nuisance values w/ no osc.
param = np.array([1.0,0.0,1.0,1.0,1.0])
#chi2 = V.LLH(param)
#print "Chi2: ",chi2


#Minimizing Chi2
nparams = 5
low_bound=np.array([0.0,-1.0,0.0,0.0,0.91])
high_bound=np.array([2.0,1.0,2.0,2.0,1.1978])
param_to_minimize=np.ones(nparams)
min_ret = V.MinLLH(param,low_bound,high_bound,param_to_minimize)
print min_ret
print "Chi2", V.GetChi2(min_ret[:-1])
"""
nuisance = min_ret[0:-1]
eprox=np.asarray( V.GetEproxEdges())
cosz=np.asarray(V.GetCosZenithEdges())

raw_expectation= V.GetPertExpectationVec(nuisance)
raw_data= V.GetDataVec()

dims=V.GetExpDims()
print "DIMS",dims

#deserialze
datasize=1
strides=[]
for i in range(0,2):
    strides.append(datasize)
    datasize*=dims[2-(i+1)]

index=0

expectation = np.zeros(dims)
expectation_nopert = np.zeros(dims)
data = np.zeros(dims)

def index(indices,sheep):
    myindex=0
    for i in range(0,2):
        myindex+=strides[i]*indices[2-(i+1)]
    #print raw_expectation[myindex] 
    return sheep[myindex]


for x in range(0,dims[0]):
    for y in range(0,dims[1]):
        myind=[x,y]
        expectation[x,y]=index(myind,raw_expectation)
        #expectation_nopert[x,y]=index(myind,raw_expectation_nopert)
        data[x,y]=index(myind,raw_data)
print min_ret[-1]

np.savez("3flav_steriledat",e=eprox,z=cosz,exp=expectation,dat=data,chi2=min_ret[-1],nuisance=nuisance)
"""





