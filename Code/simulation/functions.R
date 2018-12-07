
UBoundDF <- function(eq, params, lb.org = -Inf, ub.org = Inf){
        
        # This function calculate the lower bounds and upper bounds of 
        # the necessary condtions for that the unobserved shocks support 
        # "eq" as an equilibrium
        # 
        # Arg: (1) eq: equilibrium, a vector of 1,0; e.g eq = c(1, 0, 0, 1)
        #      (2) params: parameters in profit function
        #      (3) lb.org: original lower bound of assumed distribution
        #      (4) ub.org: orginial upper boun of assumed distribution
        #      
        # Output: a dataframe, nrow = 1
        #         "lb" & "ub" are a N column matrix, in order of firm 1, 2, 3,...,
        #         containing lower bounds and upper bounds, respectively
        #         
        
        N <- length(eq)        

        alpha <- params[grep("alpha", names(params))]
        beta.v <- params[grep("beta?", names(params))] # vector of competition effects
        
        bound <- data.frame(lb = NA, ub = NA)
        bound$lb <- matrix(lb.org, ncol = N)
        bound$ub <- matrix(ub.org, ncol = N)
        

        
        for(i in 1 : N){ 
                
                if(eq[i]) bound$ub[i] <- alpha - beta.v[-i] %*% eq[-i] 
                if(!eq[i]) bound$lb[i] <- alpha - beta.v[-i] %*% eq[-i] 
        }
        
        bound
}

LBoundDF <- function(eq, params, lb.org = -Inf, ub.org = Inf){
        
        # This function calculate the lower bounds and upper bounds of 
        # the sufficient conditions for that the unobserved shocks support "eq" as 
        # the unique equilibirum,  which is surviving iterative elimination of 
        # dominated strategy (IEDS), all possible orders of IEDS are considered
        # 
        # 
        # Arg: (1) eq: equilibrium, a vector of 1,0; e.g eq = c(1, 0, 0, 1)
        #      (2) params: parameters in profit function
        #      (3) lb.org: original lower bound of assumed distribution
        #      (4) ub.org: orginial upper boun of assumed distribution
        #      
        # Output: a dataframe, nrow = no. of permutations
        #         "lb" & "ub" are a N column matrix, in order of firm 1, 2, 3,...,
        #         containing lower bounds and upper bounds, respectively
        #         
        
        N <- length(eq)        

        gnr.bound <- matrix(rep(c(0, 1), each = NP), ncol = 2)
        
        alpha <- params[grep("alpha", names(params))]
        beta.v <- params[grep("beta?", names(params))] # vector of competition effects
        
        bound <- data.frame(orderid = 1 : NP)
        bound$order <- permutations(n = N, r = N, repeats.allowed = F)
        bound$lb <- matrix(lb.org, nrow = NP, ncol = N)
        bound$ub <- matrix(ub.org, nrow = NP, ncol = N)
        
        # matrix of entry decisions, ordered according to permutation
        entry <- matrix(eq[c(t(bound$order))], ncol = N, byrow = T)
        
        
        for(i in 1 : nrow(entry)){ 
                
                entrant <- which(entry[i, ] == 1)
                
                for(j in 1 : N){
                        
                        id <- bound$order[i, j] # firm id
                        
                        #id of entrants ordered before
                        entrant.before <- bound$order[i, entrant[entrant < j]]
                        #id of all firms ordered after
                        all.after <- bound$order[i, j : N][-1]
                        
                        #positions of firms ordered before with a larger id
                        pos.id.large <- which(bound$order[i, ] > id) 
                        id.large.before <- pos.id.large[pos.id.large < j]
                        
                        if (entry[i, j] == 1){
                                bound$ub[i, id] <- 
                                        alpha - sum(beta.v[c(entrant.before, all.after)])
                                
                                if (length(id.large.before) > 0){
                                        
                                        # latest firm with a larger id
                                        pos.before <- max(id.large.before) 
                                        bound$lb[i, id] <- 
                                                min(bound$ub[bound$order[, pos.before] == id])
                                }
                                
                                
                        } else{
                                bound$lb[i, bound$order[i, j]] <- 
                                        alpha - sum(beta.v[entrant.before])
                                
                                if (length(id.large.before) > 0){
                                        
                                        # latest firm with a larger id
                                        pos.before <- max(id.large.before) 
                                        bound$ub[i, id] <- 
                                                max(bound$lb[bound$order[, pos.before] == id])
                                }
                        }
                }
        }
        
        bound
}


ProbBounds <- function(params){
        
        # This function predicts lower bounds and upper bounds for probabilities of 
        # each possible market structure
        # 
        # Args: (1) params: vector of parameters
        # 
        # Output: dataframe with possible market structures, & their corresponding 
        # lower bounds, observed probilities, upper bounds
        

        
        prob.table <- data.frame(market.struct = matrix(numeric(0), ncol = N), 
                                 lb.prdt = numeric(0), 
                                 prob.obs = numeric(0), 
                                 up.prdt = numeric(0))
        
        for(i in 0 : 1){
                for(j in 0 : 1){
                        for(k in 0 : 1){
                                
                                if(eq == T){
                                        
                                        prob.obs <- get(paste0("P", i, j, k))
                                        
                                        UB <- UBoundDF(eq = c(i, j, k), 
                                                       params,
                                                       lb.org = 0, 
                                                       ub.org = 1)
                                        
                                        prob.ub <- apply(UB$ub - UB$lb, 1, prod)
                                        
                                        LB <- LBoundDF(eq = c(i, j, k), 
                                                       params, 
                                                       lb.org = 0, 
                                                       ub.org = 1)
                                        
                                        prob.lb <- sum(apply(LB$ub - LB$lb, 1, prod)) 
                                        
                                        prob.table <- rbind(prob.table, 
                                                            data.frame(market.struct = matrix(c(i,j,k), nrow = 1), 
                                                                       lb.prdt = prob.lb,
                                                                       prob.obs = prob.obs,
                                                                       ub.prdt = prob.ub))
                                        
                                }
                        }
                }
        }
        
        prob.table
}
