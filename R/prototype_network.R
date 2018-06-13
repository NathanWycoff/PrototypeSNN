#!/usr/bin/Rscript
#  prototype_network.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 06.03.2018

## A network of prototype neurons.
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

y1 <- 1
v1 <- 5
y2 <- 1
v2 <- 0

y1s <- rep(NA, length(ts))
v1s <- rep(NA, length(ts))
y2s <- rep(NA, length(ts))
v2s <- rep(NA, length(ts))

alpha_1 = 5
alpha_2 = 10
alpha_3 <- 0


first_reached <- NA

i <- 0
for (t in ts) {
    i <- i + 1

    y1 <- y1 + h * yp(y1, v1, t)
    v1 <- v1 + h * vp(y1, v1, t, a = 1, b = 100, 
                    alpha_1 = alpha_1, alpha_2 = alpha_2, alpha_3)
    y2 <- y2 + h * yp(y2, v2, t)
    v2 <- v2 + h * vp(y2, v2, t, a = 1, b = 100, 
                    alpha_1 = alpha_1, alpha_2 = alpha_2, alpha_3) + h*v1

    if (is.na(first_reached) && max(y1, y2) > alpha_1) {
        first_reached <- t
    }

    if (min(y1, y2) < 0) {
        break
    }

    y1s[i] <- y1
    v1s[i] <- v1
    y2s[i] <- y2
    v2s[i] <- v2
}

par(mfrow=c(1,2))
#plot(ts[1:2000], ys[1:2000])
plot(ts, y1s)
abline(h = alpha_1, col = 'green')
abline(h = alpha_2, col = 'red')
abline(h = alpha_3, col = 'blue')
abline(v = first_reached)
plot(ts, y2s)
abline(h = alpha_1, col = 'green')
abline(h = alpha_2, col = 'red')
abline(h = alpha_3, col = 'blue')
abline(v = first_reached)
