birthdeath.tree.timevar <- function (f.lamb, f.mu, lamb_par, mu_par, time.stop = 0, taxa.stop = 0, return.all.extinct=TRUE, number.species = FALSE, prune.extinct=TRUE,nmax.sp=1000000)
# start with only 1 lineage, not 2.
# allow rate variation in time as specified by f.lamb, f.mu and the parameters lamb_par and mu_par
# the parameters of the time variation are given from the past to the present
# extinct lineages can be pruned. Also return number of lineages through time.
# Note: the resulting tree has no root length. This root length is given by the first time when there is an event, i.e. $times[2]
# Written for a time.stop criteria as of now.
# Note1: do not use with the taxa.stop criterion
# with number.species not FALSE, return only the phylogenies with exactly number.species


{
    if (number.species==FALSE)
    {
    if (time.stop == 0 & taxa.stop == 0)
        stop("Must have stopping criterion\n")
            while (1) {
            	nblineages<-c(1)
            	times<-c(0)
            	b<-f.lamb(0,lamb_par)
            	d<-f.mu(0,mu_par)
            	dt <- rexp(1,(b + d))
            	#print(dt)
            	t<-dt
            	if (time.stop)
    				if (t >= time.stop) {
    					t <- time.stop
    					alive<-1
    					times<-c(times,t)
    					nblineages<-c(nblineages,1)
    					break}

            	r <- runif(1)
            	if (r>b/(b + d))
            	{
            		#print("die")
            		times<-c(times,dt)
            		nblineages<-c(nblineages,0)
            		alive<-rep(FALSE,1)}
            		else
            		{
    				edge <- rbind(c(1, 2), c(1, 3))
    				edge.length <- rep(NA, 2)
    				stem.depth <- rep(t, 2)
    				alive <- rep(TRUE, 2)
    				times<-c(times,dt)
    				nblineages<-c(nblineages,sum(alive))
    				next.node <- 4
    				repeat {
    				if (taxa.stop)
    				if (sum(alive) >= taxa.stop)
    				break
    				if (sum(alive) == 0)
    				break
    				b<-f.lamb(t,lamb_par)
        	    	d<-f.mu(t,mu_par)
    				dt <- rexp(1, sum(alive) * (b + d))
    				t <- t + dt
    				if (time.stop)
    				if (t >= time.stop) {
    					t <- time.stop
    					times<-c(times,t)
    					nblineages<-c(nblineages,sum(alive))
    					break}
    					r <- runif(1)
    					if (r <= b/(b + d)) {
    						#print("speciation")
    						random_lineage <- round(runif(1, min = 1, max = sum(alive)))
    						e <- matrix(edge[alive, ], ncol = 2)
    						parent <- e[random_lineage, 2]
    						alive[alive][random_lineage] <- FALSE
    						edge <- rbind(edge, c(parent, next.node), c(parent,next.node + 1))
    						next.node <- next.node + 2
    						alive <- c(alive, TRUE, TRUE)
    						stem.depth <- c(stem.depth, t, t)
    						#print(stem.depth)
    						x <- which(edge[, 2] == parent)
    						edge.length[x] <- t - stem.depth[x]
    						edge.length <- c(edge.length, NA, NA)
    						#print(edge.length)
    						times<-c(times,t)
    						nblineages<-c(nblineages,sum(alive))
    						if(nblineages[length(nblineages)]>=nmax.sp){print("tree is too big");break}
    						}
    						else {
    							random_lineage <- round(runif(1, min = 1, max = sum(alive)))
    							edge.length[alive][random_lineage] <- t - stem.depth[alive][random_lineage]
    							alive[alive][random_lineage] <- FALSE
    							#times<-c(times,t)
    							#nblineages<-c(nblineages,sum(alive))
    							}
    							}}
    							if (return.all.extinct == TRUE | sum(alive) > 0)
    							{
    								#print("return.tree")
    								break}
    							}}

  ######## case when the number of species is fixed #####################################################
  ######## not implemented with time-variable rates #####################################################
   	else
   	{
   		 if (time.stop == 0 & taxa.stop == 0)
        stop("Must have stopping criterion\n")
            while (1) {
            	nblineages<-c(1)
            	times<-c(0)
            	dt <- rexp(1,(b + d))
            	#print(dt)
            	t<-dt
            	#if (time.stop)
    			#	if (t >= time.stop) {
    			#		t <- time.stop
    			#		alive<-1
    			#		times<-c(times,t)
    			#		nblineages<-c(nblineages,1)
    			#		break}

            	r <- runif(1)
            	if (r>b/(b + d))
            	{
            		times<-c(times,dt)
            		nblineages<-c(nblineages,0)
            		alive<-rep(FALSE,1)}
            		else
            		{
    				edge <- rbind(c(1, 2), c(1, 3))
    				edge.length <- rep(NA, 2)
    				stem.depth <- rep(t, 2)
    				alive <- rep(TRUE, 2)
    				times<-c(times,dt)
    				nblineages<-c(nblineages,sum(alive))
    				next.node <- 4
    				repeat {
    				if (taxa.stop)
    				if (sum(alive) >= taxa.stop)
    				break
    				if (sum(alive) == 0)
    				break
    				dt <- rexp(1, sum(alive) * (b + d))
    				t <- t + dt
    				if (time.stop)
    				if (t >= time.stop) {
    					t <- time.stop
    					times<-c(times,t)
    					nblineages<-c(nblineages,sum(alive))
    					break}
    					r <- runif(1)
    					if (r <= b/(b + d)) {
    						#print("speciation")
    						random_lineage <- round(runif(1, min = 1, max = sum(alive)))
    						e <- matrix(edge[alive, ], ncol = 2)
    						parent <- e[random_lineage, 2]
    						alive[alive][random_lineage] <- FALSE
    						edge <- rbind(edge, c(parent, next.node), c(parent,next.node + 1))
    						next.node <- next.node + 2
    						alive <- c(alive, TRUE, TRUE)
    						stem.depth <- c(stem.depth, t, t)
    						#print(stem.depth)
    						x <- which(edge[, 2] == parent)
    						edge.length[x] <- t - stem.depth[x]
    						edge.length <- c(edge.length, NA, NA)
    						#print(edge.length)
    						times<-c(times,t)
    						nblineages<-c(nblineages,sum(alive))}
    						else {
    							random_lineage <- round(runif(1, min = 1, max = sum(alive)))
    							edge.length[alive][random_lineage] <- t - stem.depth[alive][random_lineage]
    							alive[alive][random_lineage] <- FALSE
    							#times<-c(times,t)
    							#nblineages<-c(nblineages,sum(alive))
    							}
    							}}
    							if (sum(alive)==number.species)
    							break
    							}
   		}

	if (sum(alive)==0) {obj<-NULL}
 	else if (sum(alive)==1) {obj<-1}
 	else
 	{
 	edge.length[alive] <- t - stem.depth[alive]
 	n <- -1
 	for (i in 1:max(edge)) {
 	if (any(edge[, 1] == i)) {
 	edge[which(edge[, 1] == i), 1] <- n
 	edge[which(edge[, 2] == i), 2] <- n
 	n <- n - 1}}
    edge[edge > 0] <- 1:sum(edge > 0)
    tip.label <- 1:sum(edge > 0)
    mode(edge) <- "character"
    mode(tip.label) <- "character"
    if(nblineages[length(nblineages)]>=nmax.sp)
    {obj<-0}
    else
    {obj <- list(edge = edge, edge.length = edge.length, tip.label = tip.label)
    class(obj) <- "phylo"
    obj <- old2new.phylo(obj)}
     if (prune.extinct)
    {obj<-prune.extinct.taxa(obj)}}
    return(list("tree"=obj,"times"=times,"nblineages"=nblineages))
}
