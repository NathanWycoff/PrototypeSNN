#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  python/analytic_test_case.py Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 07.01.2018
import numpy as np
import warnings
from scipy.optimize import minimize

class cool(object):
    pass
self = cool()

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

# npot = analytic_1neuron(4, 3)
# print(npot.fire_freq/range(20))

def exp_dec(dt, tau = 0.001):
    """
    The exponential decay postsynaptic kernel

    :param dt: Time since spike
    :param tau: PSP time scale parameter
    :return: A scalar giving the the postsynaptic current.
    """
    return np.exp(-t/tau)

I_0 = 4
V_t = 3
w = 1
tau = 0.001
t_end = 20

def fp_eq(a, b, c, init = 1, thresh = 1e-4, max_iters = 100):
    """
    Solve the equation:

    x + a x^b = c

    via fixed point iteration.

    @param init Initial value
    @param Absolute precision definite convergence
    """
    fp_iter = lambda x: c - a * (x**b)

    diff = np.inf
    x = init
    it = 0
    while (diff > thresh and it < max_iters):
        xn = fp_iter(x)
        diff = abs(xn - x)
        x = xn

    if it == max_iters:
        Warning("fp_eq did not converge")

    return x


class analytic_2neuron(object):
    """
    Simulate a two neuron system analytically.

    First neuron has constant input, second neuron is connected to the first but to no input.

    Assumes initial potential is 0 for both neurons.

    :param I_0: A scalar (usually) positive real, the constant input to the neuron.
    :param V_t: A scalar (usually) positive real, the constant threshold.
    :param tau: Time constant in exponential kernel.
    :param w: A scalar real giving the strength of the connection
    :returns: A function of time giving the value of the membrane potential.

    Example useage:

    I_0 = 1
    V_t = 1.5
    tau = 1
    w = 1
    t_end = 20

    npot = analytic_2neuron(I_0, V_t, tau, w)
    npot(t_end) 
    """
    def __init__(self, I_0, V_t, tau = 0.001, w = 1):
        self.subthresh_pot = lambda t: I_0 * (1 - np.exp(-t))
        self.w = w
        self.V_t = V_t
        self.tau = tau
        self.I_0 = I_0

    def __call__(self, t_end):
        """
        Determine firing times up to time t_end
        """
        # Firing times for presynaptic neuron
        #pres_ft = [(i+1) * self.pres_fp for i in range(int(t_end / self.pres_fp)+1)]
        # Firing times for non leaky presynaptic neuron
        fp = self.V_t / self.I_0
        pres_ft = [(i+1) * fp for i in range(int(t_end / fp)+1)]

        # Move through the sim, stopping each time the network generates a spike
        t = 0
        pot = 0
        fires1 = []
        fires2 = []
        while t < t_end:
            print(t)
            # Time of next presynaptic spike:
            Fcal = pres_ft[:int(t / fp)]
            next_pre = pres_ft[int(t / fp)]

            # See if we get a postsynaptic spike before the next presynaptic one:
            # lhs is the integral of the kernel functions
            d_const = sum([np.exp(-(t - tf) / self.tau) for tf in Fcal])
            lhs = lambda x: -tau*(sum([np.exp(-(x + t - tf)/self.tau) for tf in Fcal]) - d_const)
            rhs = (self.V_t - pot) 
            func = lambda x: lhs(x) - rhs
            obj = lambda x: np.square(func(x))

            opt_ret = minimize(obj, 0, bounds = [[0, t_end]])
            next_post = t + opt_ret['x']
            val = opt_ret['fun']

            # If it won't ever fire, or the firing time is longer than the next prespike
            if (np.sqrt(val) > 1e-2) or next_post > next_pre:
                # FFW to next presynaptic spike
                pot = pot + lhs(next_pre - t)
                t = next_pre
                fires1.append(next_pre)
            else:
                print "Postsynaptic Action Potential!"
                # FFW to the next postsynaptic spike
                pot = 0
                t = next_post
                fires2.append(next_post)

        return([fires1, fires2])
