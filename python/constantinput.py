from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
import math


mesh = IntervalMesh(128, 0, 20)

W = FunctionSpace(mesh,'CG',1)


def boundary(x):
    return x[0] < DOLFIN_EPS

def indicator(v,W,threshhold):
	ind = Function(W)
	pos = np.where(ind.vector() <= threshhold)[0]
	ind.vector()[pos] = 1.0
	print(ind.vector().array())
	return ind

def state(u):
	pos = np.where(u.vector()>=threshhold)[0]
	u.vector()[pos]=0.0
	return u


V0 = Constant("0.0")
a= Constant("1.0")
I = Constant("4.0")
v_thresh = Constant("3.0")
bc = DirichletBC(W, V0, boundary)
threshhold =3.0

v = Function(W)
u = TestFunction(W)

der = v.dx(0)

weak_form  =  der*u*dx + (a*state(v)*u - I*u)*dx

solve(weak_form == 0, v, bc)

dolfin.plot(v)
plt.show()