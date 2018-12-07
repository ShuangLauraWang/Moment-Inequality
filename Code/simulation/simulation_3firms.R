# This file uses simulation to check whether it is possible to pin down the unique eq subset as a closed-form
# Date: Apr 24, 2018
# Author: Shuang Wang
# Institution: Department of Economics, Boston University

rm(list = ls()) # clear worksapce
cat("\014") #clear the console

library(bbmle)
library(rgl)
library(alabama)
library(gtools)


# Data Generating Process -------------------------------------------------

N <- 3 #no. of firms
M <- 1000 #no. of markets

cost <- matrix(runif(N * M), ncol = N) # unobservable cost of entry
alpha <- 3/4
beta <- 1/4

data <- data.frame(mktid = seq(from = 1, to = M, by = 1),
                   cost = cost) #adding multiple columns simultaneously

#solve for eq
eq <- matrix(0, M, 2^N)

g <- 1

for(i in 0 : 1){
        for(j in 0 : 1){
                for(k in 0 : 1){
                        eq[ , g] <- (-1)^i * (alpha - beta * (2 - j - k) - cost[, 1]) >= 0 &
                                (-1)^j * (alpha - beta * (2 - i - k) - cost[, 2]) >= 0 & 
                                (-1)^k * (alpha - beta * (2 - i - j) - cost[, 3]) >= 0
                        g <- g + 1
                }
        }
}

data$eq <- sapply(apply(eq == 1, 1, which), function(x){x[sample(length(x), size = 1)]})

#different numbers indicate different: 
#eqs: 1 --> (1, 1, 1), 2 --> (1, 1, 0), 3 --> (1, 0, 1), 4 --> (1, 0, 0), 
#5 --> (0, 1, 1), 6 --> (0, 1, 0), 7 --> (0, 0, 1), 8 --> (0, 0, 0)

data$y <- cbind(data$eq <= 4, data$eq %in% c(1, 2, 5, 6), data$eq %in% c(1, 3, 5, 7))

for(i in 0 : 1){
        for(j in 0 : 1){
                for(k in 0 : 1){
                        assign(paste0("P", 1 - i, 1 - j, 1 - k),
                               sum(data$eq == c(c(i, j, k) %*% 2^((N - 1) : 0) + 1))/M)
                }
        }
}

temp <- which(order == 2, T)[,]
pos2 <- temp[order(temp[, 1]), ]

#To do: 
#(1)from bounds to possiblity;
#(2)use possiblity to construct inequalities


UB <- UBoundDF(eq = c(1, 0, 1), params, lb.org = 0, ub.org = 1)
apply(UB$ub - UB$lb, 1, prod)
LB <- LBoundDF(eq = c(1, 0, 1), params, lb.org = 0, ub.org = 1)
sum(apply(LB$ub - LB$lb, 1, prod)) 



                        
params <- setNames(c(0.75, rep(0.25, 3)),
                   c("alpha", "beta1", "beta2", "beta3"))
ProbBounds(params)




