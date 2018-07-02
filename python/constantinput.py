from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import *
import time



mesh = IntervalMesh(128, 0, 20)

W = FunctionSpace(mesh,'CG',1)


def bd_func(shift):

	def boundary(x):
		return x[0] - shift < DOLFIN_EPS
	return boundary

V0 = Constant("0.0")
a= Constant("1.0")
I = Constant("4.0")


t_mesh = np.linspace(.1,10,1000)
val_list=[]

v_current=None
t_max=20
set_log_active(False)
def cost(t,tprev,alpha,v_thresh):
	global v_current
	t=t[0]

	
	time.sleep(.5)
	mesh = IntervalMesh(128, tprev, t)
	W = FunctionSpace(mesh,'CG',1)
	bc = DirichletBC(W, V0, bd_func(tprev))
	v = Function(W)
	u = TestFunction(W)

	derv = v.dx(0)

	weak_form  =  derv*u*dx + a*v*u*dx - I*u*dx

	solve(weak_form == 0, v, bc)
	t_val = Constant(str(t))
	alpha  = Constant("0.01")
	v_current=v
	obj =  - alpha*ln(v)*dx(domain=mesh) - alpha*ln(v_thresh - v)*dx(domain=mesh)
	val = assemble(obj) - t
	# print("---")
	# print(tprev)
	# print(t)
	# print(val)
	# print("---")
	#print(val)
	if(math.isnan(val)):
		return 0
	else:
		return val
def grad(t,tprev,alpha,v_thresh):
	global v_current
	t=t[0]
	val = -1.0 - alpha*np.log(v_current(t)) - alpha*np.log(v_thresh  -v_current(t))
	#print(val)
	if(math.isnan(val)):
		return 0
	else:
		return val
def H(t,tprev,alpha,v_thresh):
	global v_current
	t=t[0]
	mesh = IntervalMesh(128, tprev, t)
	W = FunctionSpace(mesh,'CG',1)
	der = v_current.dx()
	val = -alpha*(der(t)/v_current(t)) + alpha*(der(t)/(v_thresh-v_current(t) ))
	#print(val)
	return val

t_cur= 0
alpha=.000001
v_thresh=3.0
t_hit=  []
while t_cur<=20:
	bounds =((t_cur,20),)
	x0 = np.array([t_cur+.01])
	cons = ({'type': 'ineq','fun' : lambda x: np.array([x[0] - (t_cur+.001)]),'jac' : lambda x: np.array([1.0])})
	res = minimize(cost, x0, method='SLSQP', jac=grad, bounds=bounds,constraints=cons,args=(t_cur,alpha,v_thresh))
	val = t_cur
	t_cur = res.x[0]
	if(t_cur == t_max):
		break
	t_hit.append(t_cur)


print(t_hit)

# for t in t_mesh:
# 	print(str(t))
# 	mesh = IntervalMesh(128, 0, t)

# 	W = FunctionSpace(mesh,'CG',1)
# 	bc = DirichletBC(W, V0, boundary)
# 	v = Function(W)
# 	u = TestFunction(W)

# 	derv = v.dx(0)

# 	weak_form  =  derv*u*dx + a*v*u*dx - I*u*dx

# 	solve(weak_form == 0, v, bc)
# 	t_val = Constant(str(t))
# 	alpha  = Constant("0.01")

# 	obj =  - alpha*ln(v)*dx(domain=mesh) - alpha*ln(3.0 - v)*dx(domain=mesh)
# 	val = assemble(obj) - t
# 	val_list.append(val)



# plt.plot(t_mesh,val_list)
# plt.show()

# print(val_list)
# p=Point(1.39)
# print(v(p))
# dolfin.plot(v)
# plt.show()