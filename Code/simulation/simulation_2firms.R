# This file uses simulation to check whether it is possible to estimate a Ciliberto
# & Tamer (2009) style entry model with revealed preference moment inequalities
# Date: Apr 24, 2018
# Author: Shuang Wang
# Institution: Department of Economics, Boston University

rm(list = ls()) # clear worksapce
cat("\014") #clear the console

library(bbmle)
library(rgl)
library(alabama)

# Data Generating Process -------------------------------------------------

N <- 2 #no. of firms
M <- 1000 #no. of markets

set.seed(1)
cost <- matrix(runif(N * M), ncol = N) # unobservable cost of entry
alpha <- 3/4
beta <- 1/2

data <- data.frame(mktid = seq(from = 1, to = M, by = 1),
                   cost1 = cost[, 1],
                   cost2 = cost[, 2]) #adding multiple columns simultaneously

#solve for eq
eq <- matrix(0, M, 2^2)
k <- 1


for(i in 0 : (N - 1)){
        for(j in 0 : (N - 1)){
                eq[ , k] <- (-1)^i * (alpha - beta * (1 - j) - cost[, 1]) > 0 &
                        (-1)^j * (alpha - beta * (1 - i) - cost[, 2]) > 0
                k <- k + 1
        }
}

data$eq <- sapply(apply(eq == 1, 1, which), function(x){x[sample(length(x), size = 1)]})

#different numbers indicate different eqs: 1 --> (1, 1), 2 --> (1, 0), 3 --> (0, 1), 4 --> (0, 0)


y <- matrix(0, nrow = M, ncol = N)
data$y1 <- data$eq <= 2
data$y2 <- data$eq == 1 | data$eq == 3

P11 <- sum(data$eq == 1)/M
P10 <- sum(data$eq == 2)/M
P01 <- sum(data$eq == 3)/M
P00 <- sum(data$eq == 4)/M


# Estimation --------------------------------------------------------------

obj <- function(params){
        
        alpha <- params[1]
        beta <- params[2]

        norm(min((alpha - beta) * P11 - (alpha - beta)^3/2, 0), "2") +
                norm(min(alpha * P10 -
                                 (1 - (alpha - beta)) * (alpha - beta)^2/2 -
                                 (1 - alpha) * (2 * alpha - beta) * beta/2, 0),
                     "2") +
                norm(max((alpha - beta) * P01 -
                                 alpha * (1 - alpha^2)/2 -
                                 (alpha - beta) * (2 * alpha - beta) * beta/2, 0),
                     "2") +
                norm(max(alpha * P00 - (1 - alpha) * (1 - alpha^2)/2, 0), "2") 
                # norm(min((alpha - beta) * P11 +
                #                  alpha * (1 - alpha^2)/2 -
                #                  (alpha - beta) * (2 * alpha - beta) * beta/2 +
                #                  alpha * P10 +
                #                  (1 - alpha) * (1 - alpha^2)/2, 0), "2") +
                # norm(max((alpha - beta) * P01 +
                #                  (alpha - beta)^3/2 +
                #                  alpha * P00 +
                #                  (1 - (alpha - beta)) * (alpha - beta)^2/2 -
                #                  (1 - alpha) * (2 * alpha - beta) * beta/2, 0), "2") +
                # norm(P11 - (alpha - beta)^2, "2") +
                # norm(P00 - (1 - alpha)^2, "2")

}

hin <- function(params){
        alpha <- params[1]
        beta <- params[2]
        
        h <- rep(NA, 1)
        
        h[1] <- alpha - beta
        h[2] <- 1 - alpha
        h[3] <- 1 - beta
        h
}

argnames <- c("alpha", "beta")

parnames(obj) <- argnames
parnames(hin) <- argnames

start <- setNames(c(0, 0), argnames)


auglag(fn = obj, par = start, hin = hin)
optim(fn = obj, par = start)


x <- 1 : 100/100
y <- 1 : 100/100
z <- apply(cbind(rep(x, each = 100), rep(y, times = 100)), 1, FUN = obj)


zlim <- range(z)
zlen <- round((zlim[2] - zlim[1]) * 100) + 1

colorlut <- terrain.colors(zlen)

col <- colorlut[(round(z * 100 - zlim[1]) + 1 )]

surface3d(x, y, z, color = col)
grid3d(c("x", "y+", "z"), n =20)


# graph
x <- seq(0, 1, by = 0.01)
y <- seq(0, 1, by = 0.01)
df <- setNames(expand.grid(x, y), c("x", "y"))
df <- transform(df, ueq = ((x - y) * P11 - (x - y)^3/2 >= 0) 
                & (x * P10 - (x - y)^2 * (1 - (x - y))/2 - (1 - x) * (2*x - y) * y/2 >= 0) 
                & ((x - y) * P01 - x * (1 - x^2)/2 - (x - y) * (2 * x - y) * y/2 <= 0) 
                & (x * P00 - (1 - x^2) * (1 - x)/2 <= 0)  
                & ((x - y) * P11 + 
                           x * (1 - x^2)/2 + (x - y) * (2 * x - y) * y/2 +
                           x * P10 +
                            (1 + x)/2 * P00 >= 0.5) 
                & ((x - y) * P01 + 
                           (x - y)/2 * P11 + 
                           x * P00 + 
                           (x - y)^2 * (1 - (x - y))/2 + (1 - x) * (2*x - y) * y/2 <= 0.5) 
                & (x - y >= 0)
                        )
summary(df)
df$color <- ifelse(df$ueq == TRUE, "lightblue", "orange")
with(df[df$ueq == TRUE, ], 
     plot(x = x, 
          y = y, 
          col=color, 
          type = "p",
          xlim = c(0, 1),
          ylim = c(0, 1)))



