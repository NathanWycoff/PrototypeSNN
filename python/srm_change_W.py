#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  python/srm_change_W.py Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 12.11.2018

## An arbitrary topology SRM model which allows for changing weights at each iteration.


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def W_func_sim(N, Ni, W, U, BETA, v_thresh = 3, v_reset = 0, t_eps = 1e-2, t_end = 10):
    """
    Filter-based simulation with functional-valued weight matrices, time stepping by t_eps.

    Assumes instantaneous resets (reset function is -nu*delta, where this is Dirac's delta), and a postsynaptic kernel function of exp(-dt).

    :param N: The number of neurons to simulate
    :param Ni: The number of input functions
    :param W: A function of a single variable, time, which returns a (possibly asymmetric) matrix of size (N,N), giving the connection weights of the simulated neurons at the time specified by the argument to the function.
    :param U: A function of a single variable, time, which returns a matrix of size (Ni,N) giving the connection strength from each source function to each neuron in the network at the time passed to the function.
    :param BETA: A function of a single variable, time, which returns a vector of size (Ni) giving the result of evaluating each source functional at time specified by the input.
    :param v_thresh: A real scalar: if the potential of a neuron exceeds this limit, a firing event is triggered.
    :param v_reset: A real scalar: what is the initial potential of each neuron, as well as what the potential is reset to after a firing event.
    :param t_eps: The time stepping constant.
    :param t_end: The length of the simulation. The simulation thus involves t_end / t_eps many iterations.
    """

    # Fixed postsynaptic functions for now.
    rho = lambda dt: -v_thresh if dt > 0 else 0.0
    alpha = lambda dt: (1 - np.exp(-dt)) if dt > 0 else 0.0

    t_steps = int(np.ceil(t_end/t_eps))

    # Define firing time storeage
    Fcal = [[] for _ in range(N)]

    ALPHA = np.zeros(shape=[t_steps,N])
    RHO = np.zeros(shape=[t_steps,N])
    Vs = np.empty(shape=[t_steps, N])

    t = 0
    for ti in range(t_steps):
        # Calculate the ALPHA, the postsynaptic potential, for each neuron
        for n in range(N):
            for tf in Fcal[n]:
                ALPHA[ti,n] += alpha(t - tf)

        # Calculate the RHO, the refractory contribution, for each neuron
        for n in range(N):
            for tf in Fcal[n]:
                RHO[ti,n] += rho(t - tf)

        # Sum weighted contributions from in network neurons 
        Vs[ti,:] = np.dot(W(t).T, ALPHA[ti,]) + RHO[ti,] 

        # Add on contributions from input neurons
        Vs[ti,:] += np.dot(U(t).T, BETA(t))

        # Check for firing events
        for n in range(N):
            if (Vs[ti,n] > v_thresh):
                Fcal[n].append(t)

        # Update time step
        t += t_eps

    ret = {'Fcal' : Fcal, 'Vs' : Vs}
    return ret

## Example inputs
# Define a weight matrix
np.random.seed(123)
N = 3# Number of neurons
Ni = 1# Number of input neurons
W_0 = np.random.normal(size=[N,N])
W = lambda t: W_0
U_0 = np.abs(np.random.normal(size=[Ni,N]))
U = lambda t: U_0

BETA = lambda t: np.zeros(shape=[Ni]) + t

W_func_sim(W, U, BETA)
