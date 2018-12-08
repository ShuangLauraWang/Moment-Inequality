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
library(profvis)


codewd <- "~/OneDrive - Boston University/Research/Moment Inequality/Code/simulation/"


source(paste0(codewd, "functions.R"))



# Data Generating Process -------------------------------------------------

N <- 2 #no. of firms
M <- 1000 #no. of markets

cost <- matrix(runif(N * M), ncol = N) # unobservable cost of entry
alpha <- 3/4
beta <- 1/2

data <- data.frame(mktid = seq(from = 1, to = M, by = 1),
                   cost = cost) #adding multiple columns simultaneously

#solve for eq


entry.possible <- sapply(1 : 2^N, FUN = function(x){rev(as.integer(intToBits(x - 1))[1 : N])})


eq <- apply(entry.possible, 2, FUN = function(x){
        
        temp <- sapply(1 : N, function(y){
                (-1)^x[y] *  (alpha - beta * (N - 1 - sum(x[-y])) - cost[, y]) >= 0 
        })
        
        apply(temp, 1, prod)
        
})

data$eq <- sapply(apply(eq == 1, 1, which), function(x){x[sample(length(x), size = 1)]})


# Counstruct UB and LB for the observed probs ----------------------------------


invisible(apply(entry.possible, 2, FUN = function(x){
        assign(paste0("P", paste(1 - x, collapse = "")),
               sum(data$eq == c(x %*% 2^((N - 1) : 0) + 1))/M, envir = .GlobalEnv)
}))



                        
params <- setNames(c(alpha, rep(beta, N)),
                   c("alpha", paste0("beta", 1 : N)))
ProbBounds(params)





