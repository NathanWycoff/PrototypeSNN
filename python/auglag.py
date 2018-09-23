from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import *
import time


def bd_func(shift):

	def boundary(x):
		return x[0] - shift < DOLFIN_EPS
	return boundary


val_list=[]

v_current=None
t_max=20
set_log_active(False)
def cost(t,tprev,alpha,v_thresh,mu):
	global v_current
	global mesh
	V0 = Constant("0.0")
	a= Constant("0.0")
	I = Constant("1.0")
	t1=t[0]
	s1=t[1]
	s2=t[2]
	#print(t1)
	print(mu)
	v_thresh = Constant(str(v_thresh))
	#print(tprev,t1)
	mesh = IntervalMesh(64, tprev, t1)
	W = FunctionSpace(mesh,'CG',1)
	bc = DirichletBC(W, V0, bd_func(tprev))
	v = Function(W)
	u = TestFunction(W)

	derv = v.dx(0)

	weak_form  =  derv*u*dx + a*v*u*dx - I*u*dx

	solve(weak_form == 0, v, bc)
	#t_val = Constant(str(t))
	alpha  = Constant("0.01")
	v_current=v
	obj1 =  (v_current-v_thresh)*dx(domain=mesh)
	obj2 = v_current*dx(domain=mesh)
	val = (mu/2.0)*(assemble(obj1)+s1)**2 + (mu/2.0)*(-assemble(obj2)+s2)**2 - t1
	#print("final value")
	#print(v_current(t1))
	#print("cost" + str(mu*(assemble(obj1)+s1)**2 + mu*(-assemble(obj2)+s2)**2))
	# print("---")
	# print(tprev)
	# print(t)
	# print(val)
	# print("---")
	#print(val)
	# if(math.isnan(val)):
	# 	return 0
	# else:
	# 	return val
	return val
def cont1(v_thresh):
	global v_current

def grad(t,tprev,alpha,v_thresh,mu):
	global v_current
	global mesh
	global k_counter
	V0 = Constant("0.0")
	a= Constant("0.0")
	I = Constant("1.0")
	v_thresh = Constant(str(v_thresh))
	t1=t[0]
	s1=t[1]
	s2=t[2]
	mu=mu*k_counter
	obj1 = (v_current-v_thresh)*dx(domain=mesh)
	obj2 = (v_current)*dx(domain=mesh)
	int1=assemble(obj1)
	int2=assemble(obj2)

	val1 = -1.0 + mu*(int1 + s1)*(v_current(t1)-v_thresh) + mu*(-int2+s2)*(-v_current(t1))
	val2 = mu*(int1+s1) 
	val3 = mu*(-int2+s2)
	val=(val1,val2,val3)

	#print(val)
	arr= np.array([val1,val2,val3])
	#print(arr)
	return arr
# def H(t,tprev,alpha,v_thresh):
# 	global v_current
# 	t=t[0]
# 	mesh = IntervalMesh(128, tprev, t)
# 	W = FunctionSpace(mesh,'CG',1)
# 	der = v_current.dx()
# 	val = -alpha*(der(t)/v_current(t)) + alpha*(der(t)/(v_thresh-v_current(t) ))
# 	#print(val)
# 	return val

t_cur= 0
alpha=.000001
v_thresh=1.5
t_hit=  []
counter=1
mu=10.0
shift = .000001
while t_cur<=20:
	if(t_cur+shift <20):
		bounds =((t_cur+shift,20),(0,None),(0,None))
		x0 = np.array([t_cur+shift,0.0,0.0])
		cons = ({'type': 'ineq','fun' : lambda x: np.array([x[0] - (t_cur+.001)]),'jac' : lambda x: np.array([1.0,0,0])})
		res = minimize(cost, x0, method='SLSQP', jac=grad, bounds=bounds,constraints=cons,args=(t_cur,alpha,v_thresh,mu))
		val = t_cur

		t_cur = res.x[0]
		print(res.x[0])
		#print(t_cur)
		if(t_cur == t_max):
			break
		t_hit.append(t_cur)
		#print(t_cur)
		counter+=1
		k_counter=5
	else:
		t_hit.append(20)
		break;


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