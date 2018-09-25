#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  python/analytic_subthresh.py Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 09.23.2018

## First, solve a simple linear program with NL opt
import numpy as np
from scipy.optimize import minimize

## Minimize an LP using SLSQP
obj = lambda x: sum(x)

const1 = lambda x: 2 * x[0] + x[1] - 3
const2 = lambda x: -(x[0] - 2 * x[1] - 2)
const = [{'type' : 'ineq', 'fun' : const1}, {'type' : 'ineq', 'fun' : const2}]

bound = [0, 10]
minimize(obj, [1,1], method = 'SLSQP', bounds = [bound, bound], constraints = const)

## Minimize our leaky single neuron model
I0 = 6# Input voltage
a = 1# Decay constant
nu = 3# Firing threshold
obj = lambda t: -t

c1 = -I0 / a

const1 = lambda t: -(I0/a + c1 * np.exp(-a * t) - nu)
const = [{'type' : 'ineq', 'fun' : const1}]

bound = [0, 10]
minimize(obj, [1], method = 'SLSQP', bounds = [bound], constraints = const)

## Single Active Neuron with no leak, single output spike
nu = 3# Firing threshold
w = 1# Synaptic weight connection
obj = lambda t: -t

tau = 1# Time constant in PSP
# Integrated postsynaptic potential kernel
ipsp = lambda dt: (dt > 0) * (tau - tau * np.exp(-dt / tau))

tfs = np.linspace(1, 10, 10)

const1 = lambda t: -(w * sum([ipsp(t - tf) for tf in tfs]) - nu)
const = [{'type' : 'ineq', 'fun' : const1}]

bound = [0, 10]
minimize(obj, [1], method = 'SLSQP', bounds = [bound], constraints = const)

## Single Active Neuron with no leak, multiple output spikes
nu = 3# Firing threshold
w = 1# Synaptic weight connection
t_end = 20 # Length of simulation
obj = lambda t: -t

tau = 1# Time constant in PSP
# Integrated postsynaptic potential kernel
ipsp = lambda dt: (dt > 0) * (tau - tau * np.exp(-dt / tau))

tfs = np.linspace(1, t_end, 10)

const1 = lambda t: -(w * sum([ipsp((t + current_t) - tf) for tf in tfs]) - (its) * nu)
const = [{'type' : 'ineq', 'fun' : const1}]

current_t = 0
firing_times = []
its = 1
while current_t < t_end:
    print current_t
    bound = [0, t_end - current_t]
    ret = minimize(obj, [1], method = 'SLSQP', bounds = [bound], constraints = const)
    current_t = ret['x'][0] + current_t
    firing_times.append(current_t)
    its += 1

## Two Active Neurons with no leak, multiple output spikes
nu = 3# Firing threshold
w1 = 1# Synaptic weight connections
w2 = 1
t_end = 20 # Length of simulation
obj = lambda x: -x[0]

tau = 1# Time constant in PSP
# Integrated postsynaptic potential kernel
ipsp = lambda dt: (dt > 0) * (tau - tau * np.exp(-dt / tau))

fti = np.linspace(1, t_end, 10)

thresh1 = lambda x: \
        -(w * sum([ipsp((x[1] + current_t) - tf) for tf in fti]) - (its1) * nu)
thresh2 = lambda x: \
        -(w * sum([ipsp((x[2] + current_t) - tf) for tf in ft1]) - (its2) * nu)
lin1 = lambda x: -(x[0] - x[1])
lin2 = lambda x: -(x[0] - x[2])

consts = [thresh1, thresh2, lin1, lin2]
constd = [{'type' : 'ineq', 'fun' : const} for const in consts]

# The decision variable vector x is broken down as such:
# x = [m, ta1, ta2]
current_t = 0
ft1 = []
ft2 = []
its1 = 1
its2 = 1
while current_t < t_end:
    bound = [0, t_end - current_t]

    # Do the actual optim
    ret = minimize(obj, [1,1,1], method = 'SLSQP', \
            bounds = [bound for _ in range(3)], constraints = constd)

    # Figure out who fired (1 indexed)
    x = ret['x']
    rhs = [const(x) for const in consts]
    thresh_rhs = rhs[0:2]
    who_fired = np.argmin(thresh_rhs) + 1

    # Update the happenings
    current_t = ret['x'][0] + current_t
    if (current_t < t_end):
        if (who_fired == 1):
            print "1 fired:"
            print current_t
            ft1.append(current_t)
            its1 += 1
        else:
            print "2 fired:"
            print current_t
            ft2.append(current_t)
            its2 += 1

### One input, two hidden, one output
nu = 3# Firing threshold
w01 = 1# Synaptic weight connections
w02 = 1
w13 = 1# Synaptic weight connections
w23 = 1
fire_thresh = 1e-4# Within this counts as having fired
t_end = 20 # Length of simulation
obj = lambda x: -x[0]

n_neur = 3

tau = 1# Time constant in PSP
# Integrated postsynaptic potential kernel
ipsp = lambda dt: (dt > 0) * (tau - tau * np.exp(-dt / tau))

fti = np.linspace(1, t_end, 10)

thresh1 = lambda x: \
        -(w01 * sum([ipsp((x[1] + current_t) - tf) for tf in fti]) - (its[0]) * nu)
thresh2 = lambda x: \
        -(w02 * sum([ipsp((x[2] + current_t) - tf) for tf in fti]) - (its[1]) * nu)
thresh3 = lambda x: \
        -(w13 * sum([ipsp((x[3] + current_t) - tf) for tf in fts[0]]) + 
                w23*sum([ipsp((x[3] + current_t) - tf) for tf in fts[1]]) - (its[2]) * nu)
lin1 = lambda x: -(x[0] - x[1])
lin2 = lambda x: -(x[0] - x[2])
lin3 = lambda x: -(x[0] - x[3])

consts = [thresh1, thresh2, thresh3, lin1, lin2, lin3]
constd = [{'type' : 'ineq', 'fun' : const} for const in consts]

# The decision variable vector x is broken down as such:
# x = [m, ta1, ta2, ta3]
current_t = 0
fts = [[] for _ in range(n_neur)]
its = [1 for _ in range(n_neur)]
while current_t < t_end:
    bound = [0, t_end - current_t]

    # Do the actual optim
    ret = minimize(obj, [1,1,1,1], method = 'SLSQP', \
            bounds = [bound for _ in range(n_neur+1)], constraints = constd)

    # Figure out who fired (1 indexed)
    x = ret['x']
    rhs = [const(x) for const in consts]
    thresh_rhs = rhs[0:(n_neur)]


    fired_neurs = [i for i,x in enumerate(thresh_rhs)if x < fire_thresh] 

    current_t = ret['x'][0] + current_t
    for neur in fired_neurs:
        print "%i fired!"%neur
        print current_t
        fts[neur].append(current_t)
        its[neur] += 1
