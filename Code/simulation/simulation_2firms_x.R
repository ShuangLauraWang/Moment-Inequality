# This file uses simulation to check whether it is possible to estimate a Ciliberto
# & Tamer (2009) style entry model with revealed preference moment inequalities
# Date: Apr 24, 2018
# Author: Shuang Wang
# Institution: Department of Economics, Boston University

rm(list = ls()) # clear worksapce
cat("\014") #clear the console


library(alabama)
library(expm)

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
                                  rowSums(matrix(rep(x, each = 4), nrow = 4) == 
                                                  c(1, 2, 3, 4))/
                                          length(x)})

colnames(cond.prob.df)[colnames(cond.prob.df) == "eq"] <- "cond.prob"

cond.prob.df$bin <- 1 : nrow(cond.prob.df)
data <- merge(data, cond.prob.df, by.x = c("x1.itv", "x2.itv"))

firm1 <- data[, c("mktid", "x1", "cost1", "eq", "y1", "cond.prob","x1.itv", "bin")]
firm2 <- data[, c("mktid", "x2", "cost2", "eq", "y2", "cond.prob","x2.itv", "bin")]
firm2$cond.prob <- firm2$cond.prob[, c(1, 3, 2, 4)]
colnames(firm1) <- c("mktid", "xi", "cost", "eq", "y", "cond.prob","x.itv", "bin")
colnames(firm2) <- c("mktid", "xi", "cost", "eq", "y", "cond.prob","x.itv", "bin")

firm1$firmid <- 1
firm2$firmid <- 2

data <- rbind(firm1, firm2)
data <- data[order(data$mktid), ]

temp <- matrix(data$xi, nrow = 2) 
data$xj <- c(temp[c(2, 1), ])

ineq.fn <- function(params){
        
        alpha <- params[1]
        beta <- params[2]
        
        ineq <- matrix(0, nrow(data), 6)
        
        ineq[, 1] <- (alpha * data$xi - beta) * data$cond.prob[, 1] - 
                (alpha * data$xj - beta) *  (alpha * data$xi - beta)^2/2
        
        ineq[, 2] <- alpha * data$xi * data$cond.prob[, 2] - 
                (1 - (alpha * data$xj - beta)) *  (alpha * data$xi - beta) ^2/2 - 
                (1 - alpha * data$xj) * (2 * alpha * data$xi - beta) * beta/2
        
        ineq[, 3] <- -(alpha * data$xi - beta) * data$cond.prob[, 3] + 
                alpha * data$xj * (1 - (alpha * data$xi)^2)/2 + 
                (alpha * data$xj - beta) * (2 * alpha * data$xi - beta) * beta/2

        ineq[, 4] <- -alpha * data$xi * data$cond.prob[, 4] + 
                (1 - alpha * data$xj) * (1 - (alpha * data$xi)^2)/2

        
        ineq[, 5] <- (alpha * data$xi - beta) * data$cond.prob[, 1] +
                alpha * data$xi * data$cond.prob[, 2] + 
                alpha * data$xj * (1 - (alpha * data$xi)^2)/2 + 
                (alpha * data$xj - beta) *  (2 * alpha * data$xi - beta) * beta/2 + 
                (1 + alpha * data$xi)/2 * data$cond.prob[, 4] -
                1/2

        ineq[, 6] <- -(alpha * data$xi - beta)/2 * data$cond.prob[, 1] - 
                (1 - (alpha * data$xj - beta)) *  (alpha * data$xi - beta)^2/2 - 
                (1 - alpha * data$xj) * (2 * alpha * data$xi - beta) * beta/2 -
                (alpha * data$xi - beta) * data$cond.prob[, 3] -
                alpha * data$xi * data$cond.prob[, 4] + 
                1/2
        
        ineq[, rep(1 : ncol(ineq), each = max(data$bin))] * 
                idc[, rep(1 : ncol(idc), times = 6)]
        

}

obj <- function(params){

    ineq <- ineq.fn(params)

    ineq.mean <- colMeans(ineq)
     
    ineq.sd <- apply(ineq, 2, sd)
    
    sum(pmin(ineq.mean/ineq.sd, 0)^2)
            
}

idc <- numeric(0)

for (i in 1 : max(data$bin)){
        idc <- cbind(idc, data$bin == i)
}


optim(par = c(1, 0), fn = obj)



# Confidence Sets for True Parameters -------------------------------------

## Step 1: define a grid that will contain the confidence set

hin <- function(params){
        
        ineq <- ineq.fn(params)
        
        ineq.mean <- colMeans(ineq)
        
        h <- ineq.mean + log(N * M)

}

grid.bound <- function(){
        #output: a matrix, nrow = no. of parameters, ncol = 2, 1st column are 
        #upper bounds, 2nd column are lower bounds
        
        
        D <- matrix(0, ncol = N, nrow = 2^N)
        D[seq(1, length(D), by = 2^N + 2)] <- 1
        D[seq(2, length(D), by = 2^N + 2)] <- -1
        
        bound <- apply(D, 1, 
                       FUN = function(x){
                               temp <- auglag(par = c(0, 0), 
                                              fn = function(params){x %*% params}, 
                                              hin = hin, 
                                              control.outer = list(trace = F))$par
                               }
                       )
        
        cbind(bound[seq(1, length(bound), by = 2^N + 1)], 
              bound[seq(N + 1, length(bound), by = 2^N + 1)])
        
}

grid.bound()


## Step 2: choose a point theta. With theta, we test the null hypothesis that the vecotr theta equals the true value of theta

## Step 3: evaluate the MMM test statistics at theta:

P <- 100
Alpha <- 1 : P/P
Beta <- 1 : P/P
Theta <- setNames(expand.grid(Alpha, Beta), c("alpha", "beta"))
Q <- apply(Theta, 1, FUN = obj)

## Step 4: compute the correlation matrix of the moments evaluated at theta_p

K <- 6 * max(data$bin)
R <- 1000

CS <- numeric(0)

for (i in 1 : nrow(Theta)){
        
        ineq <- ineq.fn(Theta[i, ])
        ineq.mean <- colMeans(ineq)
        ineq.sd <- apply(ineq, 2, sd)
        omega <- cor(ineq)

        rnd.xi <- matrix(rnorm(R * K), nrow = K)
        
        omega.sqrt <- chol(omega, pivot = T)
        
        temp1 <- apply(omega.sqrt %*% rnd.xi, 2, 
                       FUN = function(x){pmin(x, 0)^2})
        
        temp2 <- sqrt(N * M) * ineq.mean/ineq.sd <= sqrt(log(N * M))

        Qr <- temp2 %*% temp1
        cv <- quantile(Qr, probs = 0.95)
        
        if (Q[i] < cv){
                CS <- rbind(CS, Q[i])
        }
        
}


ineq <- ineq.fn(params)

sigma_p <- cov(ineq)
Omega <- sqrt(diag(diag(sigma))) %*% sigma %*% sqrt(diag(diag(sigma)))


## Step 5: simulate the asymptotic distribution of Q_p
## 





## 

# 2-D Graph ---------------------------------------------------------------
# 

ineq <- ineq.fn

graph.df$is <- apply(as.matrix(Theta), 1, 
                     FUN = function(x){
                             
                     })
