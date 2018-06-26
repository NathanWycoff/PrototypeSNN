from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
import math
from fenics_adjoint import *


# def boundary(x):
#     return x[0] < DOLFIN_EPS


def custom_boundary(t0):
	def boundary(x):
		return x[0]-t0 <DOLFIN_EPS
	return boundary
def cost(t,tprev,threshhold,IC):
	mesh = IntervalMesh(128, tprev, t[0])
	CG1_elem = FiniteElement("CG", mesh.ufl_cell(), 1)
	W_elem = MixedElement([CG1_elem, CG1_elem])
	W = FunctionSpace(mesh, W_elem)







t0 = 0.0
v0 = Constant("0.0")
tf=20.0


tn = t0
tnplus=None

mesh = IntervalMesh(128, 0, 20)
CG1_elem = FiniteElement("CG", mesh.ufl_cell(), 1)
W_elem = MixedElement([CG1_elem, CG1_elem])
W = FunctionSpace(mesh, W_elem)




bc = DirichletBC(W.sub(0), v0, custom_boundary(0.0))











# while tn <tf:
