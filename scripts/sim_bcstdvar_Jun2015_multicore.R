### Script for exploring parameter space and simulating trees - BCSTDVAR
### Author: Gustavo Burin
### Date: Jun 9th, 2015

# Loading libraries
library(ape)
library(geiger)
library(laser)
library(picante)
library(foreach)
library(doMC)
library(plyr)

# Loading Morlon's simulation function
source("birthdeath.tree.timevar2.R")

# Loading numbers from trees that haven't been simulated due to error
# missing <- read.table("bigtrees.txt",header=FALSE)

# Reading done trees
#done <- dir("./trees/80/complete/",pattern="_comp.txt")
#done <- gsub(pattern="_comp.txt",replacement="",x=done)
#done <- as.integer(done)
#done <- done[order(done)]
#bigtrees <- read.table("bigtrees.txt")[,1]
#new <- read.table("new_par_list.txt")

#bigtrees <- bigtrees[-match(new[,1],bigtrees)]

# Defining number of simulations
nbsim <- 2000

# Defining parameter space
lamb.min <- 0
lamb.max <- 2
# alfa.min <- 0
# alfa.max <- 1
mu.min <- 0
mu.max <- 10
beta.min <- 0
beta.max <- 1

# Setting functions to describe how both rates change through time
# Speciation decays exponentially, extinction is constant
f.lamb <- function(x,y){y}
f.mu <- function(x,y){y[1]-(y[1]*exp(-y[2]*x))}

# Function to calculate species richness in a given point in time
rich <- function(x,parms){
    a <- parms[1]
    b <- parms[2]
    c <- parms[3]
    res <- exp((a*c*x + b*(-exp(-c*x)) - b*c*x + b)/c);return(res)
}

# Creating lists to store results
par.list <- list()
par.crash <- list()
par.ext <- list()
sim.comp <- list()
sim.mol <- list()
times <- list()

# List of numbers for labelling the trees
labels <- as.list(1:nbsim)

# Uncomment if simulating only missing trees
#labels <- as.list(missing[,1])

sim <- function(x)
    {
        #print(x)
        # Creating list to store individual trees at the time of the simulation
        simi <- list()
        simi$tree <- NULL # Just to make sure the function enters the "while" part
        # The "while" loop will randomly sample values for lambda_0, alfa and mu,
        # calculate the expected theoretical values for maximum richness and final richness
        # after losing 80% of maximum richness, and if the parameters are suitable for the
        # simulations, it then calculates the time for the simulations and then runs the
        # simulation function
        while(class(simi$tree)!="phylo")
            {
                # Sampling values for mu_init and beta
                mu.par <- c(runif(1,mu.min,mu.max),runif(1,beta.min,beta.max))
                # Sampling values for lambda
                lamb.par <- mu.par[1] - runif(1,lamb.min,lamb.max)
                # print(c(lamb.par,mu.par)) # only to help debugging
                # Calculating theoretical time of decline (sp = ext)
                time.decline <- -(log((mu.par[1]-lamb.par)/mu.par[1]))/mu.par[2]
                # Calculating maximum richness at time.decline
                max.rich <- rich(time.decline,c(lamb.par,mu.par))
                # Calculating final rich, after losing 80% of max.rich
                final.rich <- 0.2*max.rich
                # Constraining trees to have between 10 an 500 expected tips as final rich
                if(final.rich>=10 & final.rich<=500)
                    {
                        # Creating a vector of expected richness through time
                        rtt <- rich(seq(0.1,100,0.1),c(lamb.par,mu.par))
                        # Finding the time at wich the tree is expected to have dropped to
                        # 20% of max.rich
                        t.final.rich <- (time.decline*10+(sum(rtt[seq(time.decline*10,1000)]>=final.rich)))/10 # Multiplied some terms and divided everything to put on the same scale
                        # Simulating complete tree (with extinct) with sampled parameters
                        simi <- birthdeath.tree.timevar(f.lamb,f.mu,lamb.par,mu.par,t.final.rich,0,return.all.extinct=TRUE,number.species=FALSE,prune.extinct=FALSE,nmax.sp=20000)
                        # Checking if the tree is valid
                        # If simi$tree == NULL, it went extinct before time ended
                        if(is.null(simi$tree))
                            {
                                par.ext <- c(par.ext,list(c(x,lamb.par,mu.par)))
                            }
                        # If simi$tree is numeric (0), the simulation reached the maximum
                        # number of tips set on the simulation function)
                        else if(class(simi$tree)=="numeric")
                            {
                                par.crash <- c(par.ext,list(c(x,lamb.par,mu.par)))
                            }
                        # If simi$tree is a valid phylogeny, then it drops extinct taxa to
                        # check if there is more than 2 extant tips. If it does not, it
                        # counts as going extinct
                        else if(class(simi$tree)=="phylo")
                            {
                                simi.mol <- drop.extinct(simi$tree)
                                if(class(simi.mol)!="phylo")
                                    {
                                        par.ext <- c(par.ext,list(c(x,lamb.par,mu.par)))
                                        simi$tree <- NULL
                                    }
                            }
                    }
                else if(final.rich<10 | final.rich>500)
                    {
                        par.fail <- c(par.fail,list(c(x,lamb.par,mu.par)))
                    }
            }
        cat(paste(x,"OK!\n\n\n",sep=" ")) # Progress indicator
        # Creating list of complete trees
        sim.comp <- c(sim.comp,list(simi))
        # Creating list of molecular trees
        sim.mol <- c(sim.mol,list(simi.mol))
        # Creating list with tree number and parameters
        write.table(data.frame(x,lamb.par,mu.par[1],mu.par[2],time.dceline,t.final.rich,final.rich),file=paste0("./pars/realized/",x,".txt"),quote=FALSE,col.names=FALSE,row.names=FALSE)
        # Exporting trees to individual text files into Newick format
        write.tree(simi.mol,file=paste("./trees/80/molecular/",x,"_mol.txt",sep=""))
        write.tree(simi$tree,file=paste("./trees/80/complete/",x,"_comp.txt",sep=""))
        # Exporting parameter lists to individual files
        write.table(matrix(unlist(par.crash),ncol=3,byrow=TRUE),file=paste0("./pars/crash/",x,"_crash.txt"),quote=FALSE,col.names=FALSE,row.names=FALSE)
        write.table(matrix(unlist(par.ext),ncol=3,byrow=TRUE),file=paste0("./pars/ext/",x,"_ext.txt"),quote=FALSE,col.names=FALSE,row.names=FALSE)
        write.table(matrix(unlist(par.fail),ncol=3,byrow=TRUE),file=paste0("./pars/fail/",x,"_fail.txt"),quote=FALSE,col.names=FALSE,row.names=FALSE)
        # Saving RData with all data
        #save.image(file="./new_sim_apr2014_multicore_richloss.RData")
        # Returning a list with complete trees, molecular trees and parameter list
        #return(list("comp.phy"=sim.comp,"mol.phy"=sim.mol,"par_list"=par.list))
    }

# Setting the number of cores to be used
registerDoMC(50)

# Running the simulations
llply(labels,sim,.parallel=TRUE)

# Saving once again the resulting object
#save(result.list,file="./new_bcstdvar_2014oct17.RData")
