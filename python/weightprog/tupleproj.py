from dolfin import *
from dolfin_adjoint import *
import matplotlib 
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
plt.ion()
import numpy as np
import math
import time
import moola

def bd_func(shift):

	def boundary(x):
		return x[0] - shift < DOLFIN_EPS
	return boundary

V0 = Constant(0.0)


## A simple bi-neural system
##
## (I)-->(weight1)>--(1)-->(weight2)>--(2)
##

# Set up function space and test functions.
mesh = IntervalMesh(1000, 0, 2*math.pi)
W = VectorFunctionSpace(mesh, 'P', 2, dim = 2)
Ws = FunctionSpace(mesh, 'P', 2)


#W = FunctionSpace(mesh, element)
bc = DirichletBC(W, [V0, V0], bd_func(0.0))
v = Function(W)
v1, v2 = split(v)
u = TestFunction(W)
u1, u2 = split(u)
#weight= interpolate(Expression('cos(x[0])-sin(x[0])', degree=2), W)
weight1 = interpolate(Expression('cos(x[0])-sin(x[0])', degree=2), Ws)
weight2 = interpolate(Expression('0.01', degree=2), Ws)

derv1 = v1.dx(0)
derv2 = v2.dx(0)

# Set up the two problems.
weak_form1  =  derv1*u1*dx -v1*u1*dx  - weight1*u1*dx - derv2*u2*dx +v2*u2*dx  + weight2*v1*u2*dx
#weak_form2  =  derv2*u2*dx -v2*u2*dx  - weight2*v1*u2*dx

Jac1 = derivative(weak_form1, v, TrialFunction(W))

solve(weak_form1==0,v, J=Jac1,bcs=bc)
plot(v1, title="state")
plot(v2, title="state")
plt.show()

sin = Expression("sin(x[0])",degree=5)


#J = .5*(v2-sin)**2*dx + 0.000000001*weight1.dx(0)**2*dx

target_time = 4.0
thresh = Constant(2.0)
const2 = Expression(str(v2([target_time])), degree = 2)
J = (const2 - thresh)**2*dx(domain=mesh)
J=assemble(J)
print(J)


rf = ReducedFunctional(J, Control(weight1))

problem = MoolaOptimizationProblem(rf)
f_moola = moola.DolfinPrimalVector(weight1)
# solver = moola.NewtonCG(problem, f_moola, options={'gtol': 1e-5,
#                                                    'maxiter': 20,
#                                                    'display': 3,
#                                                    'ncg_hesstol': 0})
solver = moola.BFGS(problem, f_moola, options={'jtol': 0,
                                               'gtol': 1e-9,
                                               'Hinit': "default",
                                               'maxiter': 100,
                                               'mem_lim': 10})




sol = solver.solve()
f_opt = sol['control'].data
plot(f_opt, title="weight1")
plt.show()

weight1.assign(f_opt)

solve(weak_form1 == 0, v, bc)

plot(v2, title="state")
xs = np.linspace(0,2*pi)
plt.plot(xs, np.sin(xs))
plt.show()
print("final objective: " + str(assemble((v1-sin)**2*dx)))
