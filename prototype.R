#!/usr/bin/Rscript
#  ../prototype.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 05.28.2018

## A prototype SNN architecture

bump_func <- function(x) {
    x <- (x - 10)*2
    if (abs(x) < 1) return(exp(-1 / (1 - x^2)))
    return(0)
}
xs <- seq(0, 15, length.out = 400)
plot(xs, sapply(xs, bump_func))

I <- function(t) {
    10
}

# I don't think this is exactly what I had in mind.
dydt <- function(y, t) {
    -y + I(t) + bump_func(y) 
}

t_max <- 5
h <- 0.001
ts <- seq(0, t_max, by = h)

ys <- rep(NA, length(ts)+1)
ys[1] <- 2

for (ti in 1:length(ts)) {
    t <- ts[ti]
    ys[ti+1] <- ys[ti] + h * dydt(ys[ti])
}

plot(ys)
