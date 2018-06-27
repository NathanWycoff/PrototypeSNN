#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  python/test_case.py Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 06.23.2018
import numpy as np
import warnings

class analytic_1neuron(object):
    """
    Simulate a single neuron analytically for constant input, constant threshold.

    Assumes initial potential is 0.

    :param I_0: A scalar (usually) positive real, the constant input to the neuron.
    :param V_t: A scalar (usually) positive real, the constant threshold.
    :returns: A function of time giving the value of the membrane potential.

    Example useage:
    npot = analytic_1neuron(4, 3)# Create a function with potential 4 and threshold 3
    npot(1) # Evaluate that function 1 second after beginning of experiment
    """
    def __init__(self, I_0, V_t):
        if (I_0 > V_t):
            self.fire_freq = np.log(I_0 / (I_0 - V_t))
        else:
            self.fire_freq = np.inf
            warnings.warn("Neuron will never fire: threshold greater than input")

        self.subthresh_pot = lambda t: I_0 * (1 - np.exp(-t))

    def __call__(self, t):
        """
        Evaluate the membrane potential at time t.
        """
        return self.subthresh_pot(t % self.fire_freq)


    



## Setup: A single neuron, weighting the input.
## input * w_1           output
## ------>-----(Neuron 1)------->
## Goal is to get the output potential to equal a certain value at the end of the time interval.

## Baseline: finite differences to integrate, finite differneces to differentiate.
## Exact arithmetic. Slow but precise.

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

### A single neuron with constant input
t_end = Decimal('2')
dt = Decimal('1e-5') #The step size for forward Euler.
int_len = int(t_end / dt)
ts = [Decimal(i) * dt for i in range(int_len)]

#Define some parameters
v_thresh = Decimal('3')#The firing threshold.
v_reset = Decimal('0')# Where we go back to after crossing v_thresh.
const_in = Decimal('4')# The constant input current.

def forward_model(weight):
    # in general, I, input current,  may depend on time, but it is constant for this example.
    I = lambda t: const_in
    Vp = lambda t, V : -V + weight*I(t)

    # Define inital voltage
    V_0 = Decimal('0')

    Vt = [Decimal('0') for _ in range(int_len)]
    Vt[0] = V_0

    F = []

    for i, t in enumerate(ts[1:]):
        Vt[i+1] = Vt[i] + dt * Vp(t, Vt[i])

        # Reset values if necessary
        if Vt[i+1] > v_thresh:
            Vt[i+1] = v_reset
            F.append(t)

    return(Vt[-1])

V_desired = Decimal('2.5')

# Define cost to get firing time as desired.
h = dt#FD step size for derivative calculation.
def cost(weight): 
    err = (V_desired - forward_model(weight))
    return err * err

fd_grad = lambda weight: (cost(weight + h) - cost(weight)) / h

init_weight = Decimal('1')

fd_grad(init_weight)

