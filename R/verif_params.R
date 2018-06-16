#!/usr/bin/Rscript
#  R/verif_params.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 06.13.2018

## Some model parameters need to be chosen, verify them here.

# Verify solution to a simple physics problem.
pred_final_vel <- function(delta, b, Vp) sqrt(Vp^2 + 2*b*delta)

pos <- 0

dist_end <- 1
h <- 1e-3

Vp <- 2
vel <- Vp
accel <- -1
poss <- c()
ts <- c()
t <- 0
while (abs(pos) < dist_end) {
    t <- t + h
    vel <- vel + accel * h
    pos <- pos + h*vel

    poss <- c(poss, pos)
    ts <- c(ts, t)
}

plot(ts, poss)


print(pred_final_vel(dist_end, accel, Vp))
print(vel)
