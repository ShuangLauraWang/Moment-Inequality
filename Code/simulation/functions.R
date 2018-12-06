
BoundDF <- function(eq, params, lb.org = -Inf, ub.org = Inf){
        
        # This function calculate the lower bounds and upper bounds of 
        # unobserved shocks for any given equilbrium for each possible order of 
        # iterative elimination of dominated strategy
        # 
        # Arg: (1) eq: equilibrium, a vector of 1,0; e.g eq = c(1, 0, 0, 1)
        #      (2) params: parameters in profit function
        #      (3) lb.org: original lower bound of assumed distribution
        #      (4) ub.org: orginial upper boun of assumed distribution
        #      
        # Output: a dataframe, nrow = no. of permutations, 
        #         "lb" & "ub" are a N column matrix, in order of firm 1, 2, 3,...,
        #         containing lower bounds and upper bounds, respectively
        #         
        
        N <- length(eq)        
        NP <- factorial(N) # no. of permutation
        
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
