from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
import math


mesh = IntervalMesh(128, 0, 20)

W = FunctionSpace(mesh,'CG',1)


def boundary(x):
    return x[0] < DOLFIN_EPS

V0 = Constant("0.0")
a= Constant("1.0")
I = Constant("4.0")
bc = DirichletBC(W, V0, boundary)

v = Function(W)
u = TestFunction(W)

der = v.dx(0)

weak_form  =  der*u*dx + a*v*u*dx - I*u*dx

solve(weak_form == 0, v, bc)

dolfin.plot(v)
plt.show()