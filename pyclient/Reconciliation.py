#! /usr/bin/python

import Verosim
import numpy as np
v=Verosim.Verosim(4)
def owl(arg):
	sheeparg=np.asarray(arg)
	print "SHEEP: ",sheeparg[0]
	print "DOG: ",sheeparg[1]
	return float(sheeparg[0]*sheeparg[1])
v.SetDeSolver(owl)
v.OscillationProbability(3,4)

