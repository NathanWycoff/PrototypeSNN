#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  python/raster_plot.py Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 12.11.2018

import matplotlib.pyplot as plt

def raster_plot(Fcal):
 """
 Given a list of firing times, make a raster plot of neural activity.

 :param Fcal: A list of lists, each inner list gives the firing times of a particular neuron, so len(Fcal) is the number of neurons, len(Fcal[0]) is the number of times the first neuron fired, etc... Further, Fcal[0][0] (should it exist) is a double giving the first time neuron 1 fired.
 """

# Calculate axis limits, a little bit paste the min and max

min_f = min([min(x) for x in Fcal])
max_f = max([max(x) for x in Fcal])
ran = max_f - min_f
min_f -= 0.1 * ran
max_f += 0.1 * ran

N = len(Fcal)
fig, axes = plt.subplots(nrows=3, ncols=1) # two axes on figure

for n in range(N):
    Fi = Fcal[n]
    for f in range(len(Fi)):
        axes[n].axvline(x=Fi[f])

for n in range(N):
    axes[n].set_xlim(left = min_f, right = max_f)

plt.suptitle("Raster Plot for %d Neuron System"%N)

plt.show()

