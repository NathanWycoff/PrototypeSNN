#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  python/analytic_test_case.py Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 07.01.2018
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




