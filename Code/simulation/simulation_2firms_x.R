# This file uses simulation to check whether it is possible to estimate a Ciliberto
# & Tamer (2009) style entry model with revealed preference moment inequalities
# Date: Apr 24, 2018
# Author: Shuang Wang
# Institution: Department of Economics, Boston University

rm(list = ls()) # clear worksapce
cat("\014") #clear the console


library(alabama)
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

for(i in 0 : (N - 1)){
    for(j in 0 : (N - 1)){
        eq[ , k] <- (-1)^i * (alpha * x[, 1] - beta * (1 - j) - cost[, 1]) > 0 &
        (-1)^j * (alpha * x[, 2] - beta * (1 - i) - cost[, 2]) > 0
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

data$x1.itv <- floor(data$x1/0.5)
data$x2.itv <- floor(data$x2/0.5)

cond.prob.df <- aggregate(eq ~ x1.itv + x2.itv, data = data,
FUN = function(x){
    rowSums(matrix(rep(x, each = 4), nrow = 4) == c(1, 2, 3, 4))/length(x)})

colnames(cond.prob.df)[colnames(cond.prob.df) == "eq"] <- "cond.prob"

data <- merge(data, cond.prob.df, by.x = c("x1.itv", "x2.itv"))

#data$cond.prob <- data$cond.prob[cbind(1 : nrow(data), data$eq)]

data <- data[order(data$mktid), ]

obj <- function(params){
    
    alpha <- params[1]
    beta <- params[2]
    
    data$ineq <- matrix(0, nrow(data), 6)
    
    data$ineq[, 1] <- (alpha * data$x1 - beta) * data$cond.prob[, 1] - 
            (alpha * data$x2 - beta) *  (alpha * data$x1 - beta)^2/2
    data$ineq[, 1] <- data$ineq[, 1] + (alpha * data$x2 - beta) * data$cond.prob[, 1] - 
            (alpha * data$x1 - beta) *  (alpha * data$x2 - beta)^2/2
    
    data$ineq[, 2] <- alpha * data$x1 * data$cond.prob[, 2] -
            (1 - (alpha * data$x2 - beta)) *  (alpha * data$x1 - beta) ^2/2 - 
                (1 - alpha * data$x2) * (2 * alpha * data$x1 - beta) * beta/2
    data$ineq[, 2] <- data$ineq[, 2] + alpha * data$x2 * data$cond.prob[, 3] -
            (1 - (alpha * data$x1 - beta)) *  (alpha * data$x2 - beta) ^2/2 - 
                (1 - alpha * data$x1) * (2 * alpha * data$x2 - beta) * beta/2
    
    data$ineq[, 3] <- -(alpha * data$x1 - beta) * data$cond.prob[, 3] + 
            alpha * data$x2 * (1 - (alpha * data$x1)^2)/2 + 
                (alpha * data$x2 - beta) * (2 * alpha * data$x1 - beta) * beta/2
    data$ineq[, 3] <- data$ineq[, 3] - (alpha * data$x2 - beta) * data$cond.prob[, 2] + 
            alpha * data$x1 * (1 - (alpha * data$x2)^2)/2 + 
                (alpha * data$x1 - beta) * (2 * alpha * data$x2 - beta) * beta/2
    
    data$ineq[, 4] <- -alpha * data$x1 * data$cond.prob[, 4] + 
            (1 - alpha * data$x2) * (1 - (alpha * data$x1)^2)/2
    data$ineq[, 4] <- data$ineq[, 4] - alpha * data$x2 * data$cond.prob[, 4] + 
            (1 - alpha * data$x1) * (1 - (alpha * data$x2)^2)/2
    
    data$ineq[, 5] <- (alpha * data$x1 - beta) * data$cond.prob[, 1] +
            alpha * data$x1 * data$cond.prob[, 2] + 
            alpha * data$x2 * (1 - (alpha * data$x1)^2)/2 + 
                (alpha * data$x2 - beta) *  (2 * alpha * data$x1 - beta) * beta/2 + 
            (1 + data$x1)/2 * data$cond.prob[, 4] -
            1/2
    data$ineq[, 5] <- data$ineq[, 5] + (alpha * data$x2 - beta) * data$cond.prob[, 1] +
            alpha * data$x2 * data$cond.prob[, 3] + 
            alpha * data$x1 * (1 - (alpha * data$x2)^2)/2 + 
                (alpha * data$x1 - beta) *  (2 * alpha * data$x2 - beta) * beta/2 + 
            (1 + data$x2)/2 * data$cond.prob[, 4] -
            1/2
    
    data$ineq[, 6] <- -(alpha * data$x1 - beta)/2 * data$cond.prob[, 1] - 
            (1 - (alpha * data$x2 - beta)) *  (alpha * data$x1 - beta) ^2/2 - 
                (1 - alpha * data$x2) * (2 * alpha * data$x1 - beta) * beta/2 -
            (alpha * data$x1 - beta) * data$cond.prob[, 3] -
            alpha * data$x1 * data$cond.prob[, 4] + 
            1/2
    data$ineq[, 6] <- data$ineq[, 6] - (alpha * data$x2 - beta)/2 * data$cond.prob[, 1] - 
            (1 - (alpha * data$x1 - beta)) *  (alpha * data$x2 - beta) ^2/2 - 
            (1 - alpha * data$x1) * (2 * alpha * data$x2 - beta) * beta/2 -
            (alpha * data$x2 - beta) * data$cond.prob[, 2] -
            alpha * data$x2 * data$cond.prob[, 4] +
            1/2
    
    data$ineq <- data$ineq/2
            
    
    ineq.mean <- aggregate(cbind(mean = ineq) ~ x1.itv + x2.itv, 
                           data = data, 
                           FUN = mean)
    
    ineq.sd <- aggregate(cbind(sd = ineq) ~ x1.itv + x2.itv, 
                           data = data, 
                           FUN = sd)
 
    
    sum(pmin(c(as.matrix(ineq.mean[, -(1 : 2)]/ineq.sd[, -(1 : 2)])), 0)^2)
    
            
}

optim(par = c(1, 0), fn = obj)

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

# graph 2-D
a <- seq(0, 1, by = 0.01)
b <- seq(0, 1, by = 0.01)
df <- setNames(expand.grid(a, b), c("a", "b"))
df <- transform(df,
ueq = ((x - y) * P11 - (x - y)^3/2 >= 0)
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





