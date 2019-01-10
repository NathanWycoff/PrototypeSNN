from dolfin import *
from dolfin_adjoint import *
import matplotlib.pyplot as plt
import numpy as np
import math
import time
import moola






def bd_func(shift):

	def boundary(x):
		return x[0] - shift < DOLFIN_EPS
	return boundary

V0 = Constant(0.0)


mesh = IntervalMesh(1000, 0, 2*math.pi)
W = FunctionSpace(mesh,'CG',1)
bc = DirichletBC(W, V0, bd_func(0.0))
v = Function(W)
u = TestFunction(W)
#weight= interpolate(Expression('cos(x[0])-sin(x[0])', degree=2), W)
weight= interpolate(Expression('0.0', degree=2), W)




derv = v.dx(0)

weak_form  =  derv*u*dx -v*u*dx  - weight*u*dx

	
Jac = derivative(weak_form, v, TrialFunction(W))
		
solve(weak_form==0,v, J=Jac,bcs=bc)
sin = Expression("sin(x[0])",degree=5)


J = .5*(v-sin)**2*dx + .000001*weight**2*dx
J=assemble(J)



rf = ReducedFunctional(J, Control(weight))

problem = MoolaOptimizationProblem(rf)
f_moola = moola.DolfinPrimalVector(weight)
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
plot(f_opt, title="weight")
plt.show()

weight.assign(f_opt)

solve(weak_form == 0, v, bc)

plot(v, title="state")
plt.show()
print("final objective: " + str(assemble((v-sin)**2*dx)))









