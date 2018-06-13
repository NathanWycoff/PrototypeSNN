#!/usr/bin/env python

# -*- coding: utf-8 -*-
#  python/neuron.py Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 06.04.2018

## A Prototype neuron architecture.
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

def yp(y, v, t):
    return v

def vp(y, v, t):
    ret = (-a)*y + \
        abs(v) * np.exp(1/y) * (v < 0) +  \
        b * (y > alpha_1 and y < alpha_2) * (v > 0) +  \
        -b *(y > alpha_2) +  \
        -b *(y > alpha_1) * (v < 0)  

    return ret

def fyt(t, y):
    return [yp(y[0], y[1], t), vp(y[0], y[1], t)]

t_int = [0, 5]
y0 = [1, 5]

alpha_1 = 5
alpha_2 = 10
alpha_3 = 0
a = 1
b = 100

sol = solve_ivp(fyt, t_int, y0, method = 'BDF')
plt.plot(sol['y'][0,:])
plt.show()
