#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  python/lp_sens.py Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 10.14.2018

## See if we can get gradients of an LP wrt A
from scipy.optimize import linprog
import numpy as np
import matplotlib.pyplot as plt


def lp(a12):
    c = [-0.5, 3]
    A = [[-3, 1], [a12, 2]]
    b = [6, 4]
    x0_bounds = (0, None)
    x1_bounds = (0, None)
    res = linprog(c, A_ub=A, b_ub=b, bounds=(x0_bounds, x1_bounds),
                  options={"disp": True})
    return(res.fun)


h = 1e-6
a12 = 1
(lp(1 + h) - lp(1)) / h


#Get grad using formula:
B = np.array([[a12, 2], [0, 1]])
bo = [4, 0]
co = [-0.5, 0]

-np.linalg.solve(B.T, co).reshape([2,1]).dot(np.linalg.solve(B.T, bo).reshape([1,2]))


## Work on getting gradients of a linearization wrt the weights

# Check our taylor series expansion
tau = 4.20
w = 0.69
tf = 0.4
f = lambda t: tau * (1.0 - np.exp(-(t - tf) / tau))

x = np.linspace(0, 100, num = 420)
plt.plot(x, [f(xi) for xi in x])
plt.show()

t_star = 3.4

f_approx = lambda t: tau * (1 - np.exp(-(t_star - tf) / tau)) + \
        np.exp(-(t_star - tf) / tau) * (t - t_star)

f_approx(3.4)
f(3.4)

# A simple SRM model
w = 0.9
nu = 3

t_opt = lambda w: tf - tau * np.log(1 - nu / (w * tau))

# Get the derivative using finite differencing and analytic stuff
h = 1e-8
anal_grad = nu * tau / (w * (nu - tau * w))
(t_opt(w + h) - t_opt(w)) / h

t_star = t_opt(w)

a = w*np.exp(-(t_star - tf) / tau)
b = nu - w * tau + w * (t_star + tau) * np.exp(-(t_star - tf) / tau)

# Match it using my dank methods
dtda = -b / np.square(a)
dadw = np.exp(-(t_star - tf) / tau)

dtdb = 1.0/a
dbdw = -tau + (t_star + tau) * np.exp(-(t_star - tf) / tau)

dtda * dadw + dtdb * dbdw
