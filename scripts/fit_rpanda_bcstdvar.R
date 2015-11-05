library(RPANDA)
library(doMC)
library(plyr)
#library(ace)
library(ape)
library(geiger)
library(laser)
library(picante)
library(paleotree)
library(phytools)

lab.80 <- dir("../trees/80/")

pars.full <- read.csv("../bcstdvar_params.csv")

registerDoMC(40)

f.lamb <- function(x,y){y}
f.mu <- function(x,y){y[1]-(y[1]*exp(-y[2]*x))}

rpanda <- function(trees,time,par)
    {
        tr <- read.tree(paste0("../trees/",time,"/",trees,"_",time,"mol.txt"))
        pars <- c(par$lambda[trees],par$mu_0[trees],par$beta[trees])
        if(time=="80")
            {
                res <- fit_bd(tr,tot_time=par$sim_time[trees],f.lamb=f.lamb,f.mu=f.mu,lamb_par=pars[1],mu_par=pars[2:3],expo.mu=FALSE,cst.lamb=TRUE,cond="stem")
            }
        if(time=="50")
            {
                res <- fit_bd(tr,tot_time=par$time50[trees],f.lamb=f.lamb,f.mu=f.mu,lamb_par=pars[1],mu_par=pars[2:3],expo.mu=FALSE,cst.lamb=TRUE,cond="stem")
            }
        if(time=="20")
            {
                res <- fit_bd(tr,tot_time=par$time20[trees],f.lamb=f.lamb,f.mu=f.mu,lamb_par=pars[1],mu_par=pars[2:3],expo.mu=FALSE,cst.lamb=TRUE,cond="stem")
            }
        if(time=="short20")
            {
                res <- fit_bd(tr,tot_time=par$time_s20[trees],f.lamb=f.lamb,f.mu=f.mu,lamb_par=pars[1],mu_par=pars[2:3],expo.mu=FALSE,cst.lamb=TRUE,cond="stem")
            }
        write.table(data.frame(res$lamb_par[1],res$mu_par[1],res$mu_par[2]),file=paste0("./results",time,"/",trees,"_fit",time,".txt"),sep=",",row.names=FALSE,col.names=FALSE,quote=FALSE)
    }

registerDoMC(40)

llply(.data=1:2000,.fun=rpanda,time="80",par=pars,.parallel=TRUE)
llply(.data=1:2000,.fun=rpanda,time="50",par=pars,.parallel=TRUE)
llply(.data=1:2000,.fun=rpanda,time="20",par=pars,.parallel=TRUE)
llply(.data=1:2000,.fun=rpanda,time="short20",par=pars,.parallel=TRUE)
