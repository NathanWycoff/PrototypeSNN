#!/usr/bin/Rscript
#  prototype2.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 06.02.2018

## A second prototype artificial neuron architecture.

# Second order as a system of two first order ODE's
yp <- function(y, v, t) v
vp <- function(y, v, t, a = 1, b = 1, alpha_1 = 1, alpha_2 = 2, alpha_3 = 0) {
    (-a)*y + #Exponential-ish decay
        abs(v) * exp(1/y) * (v < 0) + # Floor term
        b * (y > alpha_1 && y < alpha_2) * (v > 0) + #Upshoot when crossing threshold.
        -b *(y > alpha_2) + # reverses during action potential
        -b *(y > alpha_1) * (v < 0) # causes downwards movement
}

# Simulate the system using Euler's method
t_to <- 5
h <- 1e-3#Solver step size.
ts <- seq(0, t_to, by = h)

y <- 1
v <- 5

ys <- rep(NA, length(ts))
vs <- rep(NA, length(ts))

alpha_1 = 5
alpha_2 = 10
alpha_3 <- 0


first_reached <- NA

i <- 0
for (t in ts) {
    i <- i + 1

    y <- y + h * yp(y, v, t)
    v <- v + h * vp(y, v, t, a = 1, b = 100, 
                    alpha_1 = alpha_1, alpha_2 = alpha_2, alpha_3)

    if (is.na(first_reached) && y > alpha_1) {
        first_reached <- t
    }

    if (y < 0) {
        break
    }

    ys[i] <- y
    vs[i] <- v
}

par(mfrow=c(1,2))
#plot(ts[1:2000], ys[1:2000])
plot(ts, ys)
abline(h = alpha_1, col = 'green')
abline(h = alpha_2, col = 'red')
abline(h = alpha_3, col = 'blue')
abline(v = first_reached)
plot(vs)
