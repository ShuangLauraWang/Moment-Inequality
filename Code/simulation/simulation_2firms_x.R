# This file uses simulation to check whether it is possible to estimate a Ciliberto
# & Tamer (2009) style entry model with revealed preference moment inequalities
# Date: Apr 24, 2018
# Author: Shuang Wang
# Institution: Department of Economics, Boston University

rm(list = ls()) # clear worksapce
cat("\014") #clear the console
# Data Generating Process -------------------------------------------------

N <- 2 #no. of firms
M <- 1000 #no. of markets

set.seed(1)
cost <- matrix(runif(N * M), ncol = N) # unobservable cost of entry
x <-  matrix(runif(N * M, min = 0, max = 2), ncol = N)

alpha <- 3/4
beta <- 1/2

data <- data.frame(mktid = seq(from = 1, to = M, by = 1),
                   x1 = x[, 1],
                   x2 = x[, 2],
                   cost1 = cost[, 1],
                   cost2 = cost[, 2]) #adding multiple columns simultaneously

#solve for eq
eq <- matrix(0, M, 2^N)
k <- 1

for(j in 0 : (N - 1)){
        for(i in 0 : (N - 1)){
                eq[ , k] <- (-1)^i * (alpha * x[, 1] - beta * j - cost[, 1]) < 0 & 
                        (-1)^j * (alpha * x[, 2] - beta * i - cost[, 2]) < 0
                k <- k + 1
        }
}

data$eq <- sapply(apply(eq == 1, 1, which), function(x){x[sample(length(x), size = 1)]})

#different numbers indicate different eqs: 1 --> (0, 0), 2 --> (1, 0), 3 --> (0, 1), 4 --> (1, 1)

y <- matrix(0, nrow = M, ncol = N)
data$y1 <- data$eq == 2 | data$eq == 4
data$y2 <- data$eq >= 3

P11 <- sum(data$eq == 1)/M
P10 <- sum(data$eq == 2)/M
P01 <- sum(data$eq == 3)/M
P00 <- sum(data$eq == 4)/M

data$x1.itv <- floor(data$x1/0.5)
data$x2.itv <- floor(data$x2/0.5)

cond.prob.df <- aggregate(eq ~ x1.itv + x2.itv, data = data, 
                          FUN = function(x){
                                  rowSums(matrix(rep(x, each = 4), nrow = 4) == c(1, 2, 3, 4))/length(x)})

colnames(cond.prob.df)[colnames(cond.prob.df) == "eq"] <- "cond.prob"

cond.prob.df <- merge(data, cond.prob.df, by.x = c("x1.itv", "x2.itv"))

cond.prob.df <- cond.prob.df[order(cond.prob.df$mktid), ]


data$cond.prob <- cond.prob.df$cond.prob[cbind(1 : nrow(data), data$eq)]


obj <- function(params){
        
                alpha <- params[1]
                beta <- params[2]
        
                Q <- 0
                
                sub.data1 <- data[data$eq == 1, ]
                Q <- Q + sum(pmax(alpha * sub.data1$x1 * sub.data1$cond.prob - 
                                          (1 - alpha * sub.data1$x2) * (1 - (alpha * sub.data1$x1))/2, 
                                  0)^2)
                Q <- Q + sum(pmax(alpha * sub.data1$x2 * sub.data1$cond.prob - 
                                          (1 - alpha * sub.data1$x1) * (1 - (alpha * sub.data1$x2))/2, 
                                  0)^2)
                
                sub.data2 <- data[data$eq == 2, ]
                
                Q <- Q + sum(pmin(alpha * sub.data2$x1 * sub.data2$cond.prob - 
                                          (1 - (alpha * sub.data2$x2 - beta)) * (alpha * sub.data2$x1 - beta)^2/2 - 
                                          (1 - alpha * sub.data2$x2) * ((alpha * sub.data2$x1)^2 - 
                                                                                (alpha * sub.data2$x1 - beta)^2)/2, 
                                  0)^2)
                
                
                Q <- Q + sum(pmax((alpha * sub.data2$x2 - beta) * sub.data2$cond.prob -
                                          alpha * sub.data2$x1 * (1 - (alpha * sub.data2$x2)^2)/2 -
                                          (alpha * sub.data2$x1 - beta) * ((alpha * sub.data2$x2)^2 - 
                                                                                (alpha * sub.data2$x2 - beta)^2)/2,
                                   0)^2)
               
                sub.data3 <- data[data$eq == 3, ]
                Q <- Q + sum(pmax((alpha * sub.data3$x1 - beta) * sub.data3$cond.prob -
                                          alpha * sub.data3$x2 * (1 - (alpha * sub.data3$x1)^2)/2 -
                                          (alpha * sub.data3$x2 - beta) * ((alpha * sub.data3$x1)^2 - 
                                                                                (alpha * sub.data3$x1 - beta)^2)/2,
                                  0)^2)
                
                Q <- Q + sum(pmin(alpha * sub.data3$x2 * sub.data3$cond.prob - 
                                          (1 - (alpha * sub.data3$x1 - beta)) * (alpha * sub.data3$x2 - beta)^2/2 - 
                                          (1 - alpha * sub.data3$x1) * ((alpha * sub.data3$x2)^2 - 
                                                                                (alpha * sub.data3$x2 - beta)^2)/2, 
                                  0)^2)
                
                
                sub.data4 <- data[data$eq == 4, ]
                Q <- Q + sum(pmin((alpha * sub.data4$x1 - beta) * sub.data4$cond.prob - 
                                      (alpha * sub.data4$x2 - beta) * (alpha * sub.data4$x1 - beta)^2/2, 
                                  0)^2)
                Q <- Q + sum(pmin((alpha * sub.data4$x2 - beta) * sub.data4$cond.prob - 
                                          (alpha * sub.data4$x1 - beta) * (alpha * sub.data4$x2 - beta)^2/2, 
                                  0)^2)
                return(Q/M)
}

obj.eq <- function(params){
        alpha <- params[1]
        beta <- params[2]
        
        Q <- 0
        
        sub.data1 <- data[data$eq == 1, ]
        Q <- Q + 
                sum(((1 - punif(alpha * sub.data1$x1)) * ((1 - punif(alpha * sub.data1$x2))) - 
                             sub.data1$cond.prob)^2)
        
        sub.data4 <- data[data$eq == 4, ]
        Q <- Q + sum((punif(alpha * sub.data4$x1 - beta) * punif(alpha * sub.data4$x2 - beta) - 
                              sub.data4$cond.prob)^2)
        
}

obj.eq.in <- function(params){
        
        alpha <- params[1]
        beta <- params[2]
        
        
        Q <- 0
        
        sub.data1 <- data[data$eq == 1, ]
        Q <- Q + 
                sum(((1 - punif(alpha * sub.data1$x1)) * ((1 - punif(alpha * sub.data1$x2))) - 
                              sub.data1$cond.prob)^2)

        sub.data2 <- data[data$eq == 2, ]
        
        Q <- Q + sum(pmin(alpha * sub.data2$x1 * sub.data2$cond.prob - 
                                  (1 - (alpha * sub.data2$x2 - beta)) * (alpha * sub.data2$x1 - beta)^2/2 - 
                                  (1 - alpha * sub.data2$x2) * ((alpha * sub.data2$x1)^2 - 
                                                                        (alpha * sub.data2$x1 - beta)^2)/2, 
                          0)^2)
        
        
        Q <- Q + sum(pmax((alpha * sub.data2$x2 - beta) * sub.data2$cond.prob -
                                  alpha * sub.data2$x1 * (1 - (alpha * sub.data2$x2)^2)/2 -
                                  (alpha * sub.data2$x1 - beta) * ((alpha * sub.data2$x2)^2 - 
                                                                           (alpha * sub.data2$x2 - beta)^2)/2,
                          0)^2)
        
        sub.data3 <- data[data$eq == 3, ]
        Q <- Q + sum(pmax((alpha * sub.data3$x1 - beta) * sub.data3$cond.prob -
                                  alpha * sub.data3$x2 * (1 - (alpha * sub.data3$x1)^2)/2 -
                                  (alpha * sub.data3$x2 - beta) * ((alpha * sub.data3$x1)^2 - 
                                                                           (alpha * sub.data3$x1 - beta)^2)/2,
                          0)^2)
        
        Q <- Q + sum(pmin(alpha * sub.data3$x2 * sub.data3$cond.prob - 
                                  (1 - (alpha * sub.data3$x1 - beta)) * (alpha * sub.data3$x2 - beta)^2/2 - 
                                  (1 - alpha * sub.data3$x1) * ((alpha * sub.data3$x2)^2 - 
                                                                        (alpha * sub.data3$x2 - beta)^2)/2, 
                          0)^2)
        
        sub.data4 <- data[data$eq == 4, ]
        Q <- Q + sum((punif(alpha * sub.data4$x1 - beta) * punif(alpha * sub.data4$x2 - beta) - 
                     sub.data4$cond.prob)^2)
}

optim(par = c(0, 0), fn = obj.eq)

#graph
x <- 1 : 100/100
y <- 1 : 100/100
z <- apply(cbind(rep(x, each = 100), rep(y, times = 100)), 1, FUN = obj)

df <- data.frame(x = rep(x, each = 100), y = rep(y, times = 100), z = z)


zlim <- range(z)
zlen <- round((zlim[2] - zlim[1]) * 100) + 1

colorlut <- terrain.colors(zlen)

col <- colorlut[(round(z * 100 - zlim[1]) + 1 )]

surface3d(x, y, z, color = col)
grid3d(c("x", "y+", "z"), n =20)



