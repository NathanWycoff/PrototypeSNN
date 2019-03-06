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


class Delta(UserExpression):
    def __init__(self, **kwargs):
        UserExpression.__init__(self, **kwargs)
    def eval(self, values, x):
        values[0] = my_expr

    def value_shape(self): return (1, )


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
tmax = 2*math.pi
dim = 2
mesh = IntervalMesh(1000, 0, tmax)
W = VectorFunctionSpace(mesh, 'P', 2, dim = dim)
Ws = FunctionSpace(mesh, 'P', 2)

#boundary_parts = \
#    MeshFunction('size_t', mesh, mesh.topology().dim() - 1)
#
#right = CompiledSubDomain("near(x[0], side)", side = tmax)
#right.mark(boundary_parts, 0)
#dPP = dP(subdomain_data=boundary_parts)

#W = FunctionSpace(mesh, element)
bc = DirichletBC(W, [V0 for _ in range(dim)], bd_func(0.0))
v = Function(W)
vs = split(v)
u = TestFunction(W)
us = split(u)
#weight= interpolate(Expression('cos(x[0])-sin(x[0])', degree=2), W)
weight1 = interpolate(Expression('cos(x[0])-sin(x[0])', degree=2), Ws)
weight2 = interpolate(Expression('0.01', degree=2), Ws)
assign(v1, interpolate(Expression('1.0', degree=2), Ws))
v2 = interpolate(Expression('1.0', degree=2), Ws)
v1, v2 = vs
u1, u2 = us

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

delta=Delta()

my_expr = Expression("0.0001/3.14/(((x[0])-4)^2 + 0.0001^2)", degree = 2)
my_expr = Expression("1E-2 / pi * 1.0/(pow(x[0]-4.0,2) + 1E-4)", degree = 10)

target_time = 4.0#TODO: Incorporate Target Time.
thresh = 1.0
#assign(v, [Constant(1.0), Constant(1.0)])
#assign(v1, Constant(1.0))
J = (v1 - thresh)**2*my_expr*dx + 100*weight1.dx(0)**2*dx 
J=assemble(J)
print(J)

rf = ReducedFunctional(J, Control(weight1))

problem = MoolaOptimizationProblem(rf)
f_moola = moola.DolfinPrimalVector(weight1)
solver = moola.NewtonCG(problem, f_moola, options={'gtol': 1e-5,
                                                    'maxiter': 20,
                                                    'display': 3,
                                                    'ncg_hesstol': 0})
#solver = moola.BFGS(problem, f_moola, options={'jtol': 0,
#                                               'gtol': 1e-9,
#                                               'Hinit': "default",
#                                               'maxiter': 100,
#                                               'mem_lim': 10})




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
