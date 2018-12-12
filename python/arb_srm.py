#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  python/arb_srm.py Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 12.10.2018

## An implementation of an SRM time stepper with arbitrary topology.

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Define some parameters
v_thresh = 3#The firing threshold.
v_reset = 0# Where we go back to after crossing v_thresh.
const_in = 4# The constant input current.

# in general, I, input current,  may depend on time, but it is constant for this example.
TAU = 1
alpha = lambda dt: TAU * (1 - np.exp(-dt / TAU)) if dt > 0 else 0.0
rho = lambda dt: -v_thresh if dt > 0 else 0.0

# Define inital voltage
V_0 = 0

t_eps = 1e-2
t_end = 10
t_steps = int(np.ceil(t_end/t_eps))

# Define a weight matrix
N = 3# Number of neurons
Ni = 1# Number of input neurons
W = np.random.normal(size=[N,N])
U = np.random.normal(size=[Ni,N])

# Define firing time storeage
Fcal = [[] for _ in range(N)]

# Define input firing times
n_fires = 100
Fi = [np.random.gamma(1,1,size=[n_fires])]

ALPHA = np.zeros(shape=[t_steps,N])
BETA = np.zeros(shape=[t_steps,Ni])
RHO = np.zeros(shape=[t_steps,N])

t = 0
for ti in range(t_steps):
    # Calculate the ALPHA, the postsynaptic potential, for each neuron
    for n in range(N):
        for tf in Fcal[n]:
            ALPHA[ti,n] += alpha(t - tf)

    # Calculate the BETA, the input potential
    for n in range(Ni):
        for tf in Fi[n]:
            BETA[ti,n] += alpha(t - tf)

    # Calculate the RHO, the refractory contribution, for each neuron
    for n in range(N):
        for tf in Fcal[n]:
            RHO[ti,n] += rho(t - tf)

    # Sum weighted contributions from in network neurons 
    V = np.dot(W.T, ALPHA[ti,]) + RHO[ti,] 

    # Add on contributions from input neurons
    V += np.dot(U.T, BETA[ti,])

    # Check for firing events
    for n in range(N):
        if (V[n] > v_thresh):
            Fcal[n].append(t)

    # Update time step
    t += t_eps
