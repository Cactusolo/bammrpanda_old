data <- read.csv("~/Dropbox/bamm_panda/bvardcst/results/fulldata_bvardcst.csv")

tree <- read.tree("~/Documents/temp_bammpanda/bvardcst/complete/80/6_80comp.txt")
tree.mol <- read.tree("~/Dropbox/bamm_panda/bvardcst/trees/80/6_mol80_8.162_0.921_0.128_1.164_21_.txt")

source("~/Dropbox/collaborations/Le_Codes/branching.times.with.extinction.R")


pitoco <- data$time.80[6]-max(branching.times.with.extinct(tree))


# Function to calculate species richness in a given point in time
rich <- function(x,parms){
    a <- parms[1]
    b <- parms[2]
    c <- parms[3]
    res <- exp((a*(-exp(-b*x))+a-(b*c*x))/b);return(res)
}

## trees80 <- lapply(paste0("~/Dropbox/bamm_panda/bcstdvar/trees/80/",1:2000,"_80mol.txt"),read.tree)
## trees50 <- lapply(paste0("~/Dropbox/bamm_panda/bcstdvar/trees/50/",1:2000,"_50mol.txt"),read.tree)
## trees20 <- lapply(paste0("~/Dropbox/bamm_panda/bcstdvar/trees/20/",1:2000,"_20mol.txt"),read.tree)
## trees.s20 <- lapply(paste0("~/Dropbox/bamm_panda/bcstdvar/trees/short20/",1:2000,"_short20mol.txt"),read.tree)

## cga.80 <- unlist(lapply(trees80,function(x) max(branching.times(x))))
## cga.50 <- unlist(lapply(trees50,function(x) max(branching.times(x))))
## cga.20 <- unlist(lapply(trees20,function(x) max(branching.times(x))))
## cga.s20 <- unlist(lapply(trees.s20,function(x) max(branching.times(x))))

## bcstdvar$cga.80 <- cga.80
## bcstdvar$cga.50 <- cga.50
## bcstdvar$cga.20 <- cga.20
## bcstdvar$cga.s20 <- cga.s20

source("~/Dropbox/collaborations/Le_Codes/gammaCorrected.R")

gamma.80 <- unlist(lapply(trees80,gammaCorrected))
gamma.50 <- unlist(lapply(trees50,gammaCorrected))
gamma.20 <- unlist(lapply(trees20,gammaCorrected))
gamma.s20 <- unlist(lapply(trees.s20,gammaCorrected))

bcstdvar$gamma.80 <- gamma.80
bcstdvar$gamma.50 <- gamma.50
bcstdvar$gamma.20 <- gamma.20
bcstdvar$gamma.s20 <- gamma.s20


time.decline.80 <- bcstdvar$time80 + (log((bcstdvar$mu_0.sim-bcstdvar$lambda.sim)/bcstdvar$mu_0.sim))/bcstdvar$beta.sim
time.decline.50 <- bcstdvar$time50 + (log((bcstdvar$mu_0.sim-bcstdvar$lambda.sim)/bcstdvar$mu_0.sim))/bcstdvar$beta.sim
time.decline.20 <- bcstdvar$time20 + (log((bcstdvar$mu_0.sim-bcstdvar$lambda.sim)/bcstdvar$mu_0.sim))/bcstdvar$beta.sim


bcstdvar$decline.length.80 <- time.decline.80
bcstdvar$decline.length.50 <- time.decline.50
bcstdvar$decline.length.20 <- time.decline.20

write.table(bcstdvar,file="~/Dropbox/bamm_panda/bcstdvar/results/fulldata_bcstdvar.csv",quote=FALSE,row.names=FALSE,sep=",")

par(mfrow=c(2,1))


pdf("phylo_timetravel_2.pdf")
plot(tree,show.tip.label=FALSE);axisPhylo()
abline(v=max(branching.times.with.extinct(tree)),col=2, lty=2, lwd=2)
abline(v=data$time.50[6], col=3, lty=2, lwd=2)
abline(v=data$time.20[6], col=4, lty=2, lwd=2)
abline(v=data$time.s20[6], col=6, lty=2, lwd=2)
dev.off()

plot(y=rich(x=seq(pitoco,data$time.80[6],0.01)-pitoco,parms=c(data$lambda.init.sim.80[6],data$alpha.sim.80[6],data$mu.sim.80[6])),x=seq(pitoco,data$time.80[6],0.01)-pitoco,type="l",ylab="Species Richness",xlab="Time")
abline(v=max(branching.times.with.extinct(tree)),col=2, lty=2)
abline(v=data$time.50[6], col=3, lty=2)
abline(v=data$time.20[6], col=4, lty=2)
abline(v=data$time.s20[6], col=6, lty=2)


rtt <- data.frame("rich"=rich(x=seq(pitoco,data$time.80[6],0.01)-pitoco,parms=c(data$lambda.init.sim.80[6],data$alpha.sim.80[6],data$mu.sim.80[6])),"time"=seq(pitoco,data$time.80[6],0.01)-pitoco)

library(ggplot2)
library(devtools)
source("http://bioconductor.org/biocLite.R")
biocLite("ggtree")
library(ggtree)

g1 <- ggplot(data=rtt,aes(x=time,y=rich)) + geom_line() + theme_bw() + geom_vline(xintercept=max(branching.times.with.extinct(tree)),linetype="longdash",colour="red") + geom_vline(xintercept=data$time.50[6],linetype="longdash",colour="green3") + geom_vline(xintercept=data$time.20[6],linetype="longdash",colour="blue") + geom_vline(xintercept=data$time.s20[6],linetype="longdash", colour="magenta") + annotate("text",x=max(branching.times.with.extinct(tree)) + 0.2,y=200,label="80%",colour="red") + annotate("text",x=data$time.50[6] - 0.2,y=200,label="50%",colour="green3") + annotate("text",x=data$time.20[6] - 0.2,y=200,label="20%",colour="blue") + annotate("text",x=data$time.s20[6] - 0.35,y=200,label="short20%",colour="magenta") + theme(axis.text.x=element_blank()) + ylab("Species richness") + xlab("Time")

ggsave("diversity_timetravel_2.pdf",g1)



ggsave("ntip80.pdf",ggplot(data=data,aes(x=ntip.80)) + geom_histogram(binwidth=1) + theme_bw() + xlab("Number of Species") + ylab("Count"))
ggsave("simtime.pdf",ggplot(data=data,aes(x=time.80)) + geom_histogram(binwidth=1) + theme_bw() + xlab("Tree Length (Myrs)") + ylab("Count"))
ggsave("gamma80.pdf",ggplot(data=data,aes(x=gamma.80)) + geom_histogram(binwidth=0.01) + theme_bw() + geom_vline(xintercept=0,linetype="longdash",colour="red") + xlab("Corrected Gamma Statistic") + ylab("Count"))


bvardcst <- data
bcstdvar <- read.csv("~/Dropbox/bamm_panda/bcstdvar/results/fulldata_bcstdvar.csv")


## BVARDCST mu by lambda init
bvdc <- ggplot(data=bvardcst,aes(x=lambda.init.sim.80,y=mu.sim.80)) + geom_point(aes(color=alpha.sim.80)) + theme_bw() + theme(legend.position="bottom") + scale_colour_gradient(low="#24C6DC",high="#514A9D",name="Alpha") + ylim(0,10) + xlim(0,10) + xlab("Initial Lambda") + ylab("Mu")
## BCSTDVAR mu final by lambda
bcdv <- ggplot(data=bcstdvar,aes(x=lambda.sim,y=mu_0.sim)) + geom_point(aes(color=beta.sim)) + theme_bw() + theme(legend.position="bottom") + scale_colour_gradient(low="#FFC837",high="#FF8008",name="Beta") + ylim(0,10) + xlim(0,10) + xlab("Lambda") + ylab("Final Mu")

ggsave("bvardcst_parspace.pdf",bvdc)
ggsave("bcstdvar_parspace.pdf",bcdv)





bvdc.col <- colorRampPalette(c("#24C6DC","#514A9D"))
bcdv.col <- colorRampPalette(c("#FFC837","#FF8008"))




## BVARDCST lambda final
b <- ggplot(data=bvardcst,aes(x=lambda.final.sim.80,y=mu.sim.80)) +
    geom_point(aes(colour=alpha.sim.80)) +
        theme_bw() +
            theme(legend.position="bottom") +
                scale_colour_gradient(low="#24C6DC",high="#514A9D") +
                    ylim(0,2.5) +
                        xlim(0,2.5)



## BVARDCST R initial
a <- ggplot(data=bvardcst,aes(x=lambda.init.sim.80 - mu.sim.80)) +
    geom_histogram(aes(y=..density..),fill=bvdc.col(10)[1],alpha=0.75) +
        geom_density(colour=bvdc.col(10)[10],size=1.5,fill=bvdc.col(10)[10],alpha=0.25) +
            theme_bw() +
                xlim(0,10) +
                    ylim(0,0.25) +
                        xlab("Initial Net Diversification") +
                            ylab("Probability Density")
ggsave("bvardcst_initial_r.pdf",a)

## BCSTDVAR R initial
b <- ggplot(data=bcstdvar,aes(x=lambda.sim)) +
    geom_histogram(aes(y=..density..),fill=bcdv.col(10)[1],alpha=0.75) +
        geom_density(colour=bcdv.col(10)[10],size=1.5,fill=bcdv.col(10)[10],alpha=0.25) +
            theme_bw() +
                xlim(0,10) +
                    ylim(0,0.25) +
                        xlab("Initial Net Diversification") +
                            ylab("Probability Density")
ggsave("bcstdvar_initial_r.pdf",b)

## BVARDCST R final
d <- ggplot(data=bvardcst,aes(x=lambda.final.sim.80 - mu.sim.80)) +
    geom_histogram(aes(y=..density..),fill=bvdc.col(10)[1],alpha=0.75) +
        geom_density(colour=bvdc.col(10)[10],size=1.5,fill=bvdc.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Final Net Diversification") +
                            ylab("Probability Density")
ggsave("bvardcst_final_r.pdf",d)

## BCSTDVAR R final
e <- ggplot(data=bcstdvar,aes(x=lambda.sim - mu_0.sim)) +
    geom_histogram(aes(y=..density..),fill=bcdv.col(10)[1],alpha=0.75) +
        geom_density(colour=bcdv.col(10)[10],size=1.5,fill=bcdv.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Initial Net Diversification") +
                            ylab("Probability Density")
ggsave("bcstdvar_final_r.pdf",e)

## BVARDCST Crown Group Age
g <- ggplot(bvardcst,aes(x=cga.80)) +
    geom_histogram(aes(y=..density..),fill=bcdv.col(10)[1],alpha=0.75) +
        geom_density(colour=bcdv.col(10)[10],size=1.5,fill=bcdv.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Crown Group Age") +
                            ylab("Probability Density")
ggsave("bcstdvar_final_r.pdf",g)

## BCSTDVAR Crown Group Age
h <- ggplot(data=bcstdvar,aes(x=cga.80)) +
    geom_histogram(aes(y=..density..),fill=bvdc.col(10)[1],alpha=0.75) +
        geom_density(colour=bvdc.col(10)[10],size=1.5,fill=bvdc.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Crown Group Age") +
                            ylab("Probability Density")
ggsave("bcstdvar_final_r.pdf",h)

## BVARDCST ntip 80
h <- ggplot(data=bvardcst,aes(x=ntip.80)) +
    geom_histogram(aes(y=..density..),fill=bvdc.col(10)[1],alpha=0.75,binwidth=max(bvardcst$ntip.80)/100) +
        geom_density(colour=bvdc.col(10)[10],size=1.5,fill=bvdc.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Number of Tips") +
                            ylab("Probability Density")
ggsave("bvardcst_ntip_80.pdf",h)

## BVARDCST ntip 50
h <- ggplot(data=bvardcst,aes(x=ntip.50)) +
    geom_histogram(aes(y=..density..),fill=bvdc.col(10)[1],alpha=0.75,binwidth=max(bvardcst$ntip.50)/100) +
        geom_density(colour=bvdc.col(10)[10],size=1.5,fill=bvdc.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Number of Tips") +
                            ylab("Probability Density")
ggsave("bvardcst_ntip_50.pdf",h)

## BVARDCST ntip 20
h <- ggplot(data=bvardcst,aes(x=ntip.20)) +
    geom_histogram(aes(y=..density..),fill=bvdc.col(10)[1],alpha=0.75,binwidth=max(bvardcst$ntip.20)/100) +
        geom_density(colour=bvdc.col(10)[10],size=1.5,fill=bvdc.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Number of Tips") +
                            ylab("Probability Density")
ggsave("bvardcst_ntip_20.pdf",h)

## BVARDCST ntip short20
h <- ggplot(data=bvardcst,aes(x=ntip.s20)) +
    geom_histogram(aes(y=..density..),fill=bvdc.col(10)[1],alpha=0.75,binwidth=max(bvardcst$ntip.s20)/100) +
        geom_density(colour=bvdc.col(10)[10],size=1.5,fill=bvdc.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Number of Tips") +
                            ylab("Probability Density")
ggsave("bvardcst_ntip_s20.pdf",h)



## BVARDCST gamma 80
h <- ggplot(data=bvardcst,aes(x=gamma.80)) +
    geom_histogram(aes(y=..density..),fill=bvdc.col(10)[1],alpha=0.75) +
        geom_density(colour=bvdc.col(10)[10],size=1.5,fill=bvdc.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Corrected Gamma") +
                            ylab("Probability Density")
ggsave("bvardcst_gamma_80.pdf",h)

## BVARDCST gamma 80
h <- ggplot(data=bvardcst,aes(x=gamma.50)) +
    geom_histogram(aes(y=..density..),fill=bvdc.col(10)[1],alpha=0.75) +
        geom_density(colour=bvdc.col(10)[10],size=1.5,fill=bvdc.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Corrected Gamma") +
                            ylab("Probability Density")
ggsave("bvardcst_gamma_50.pdf",h)

## BVARDCST gamma 20
h <- ggplot(data=bvardcst,aes(x=gamma.20)) +
    geom_histogram(aes(y=..density..),fill=bvdc.col(10)[1],alpha=0.75) +
        geom_density(colour=bvdc.col(10)[10],size=1.5,fill=bvdc.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Corrected Gamma") +
                            ylab("Probability Density")
ggsave("bvardcst_gamma_20.pdf",h)

## BVARDCST gamma short20
h <- ggplot(data=bvardcst,aes(x=gamma.s20)) +
    geom_histogram(aes(y=..density..),fill=bvdc.col(10)[1],alpha=0.75) +
        geom_density(colour=bvdc.col(10)[10],size=1.5,fill=bvdc.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Corrected Gamma") +
                            ylab("Probability Density")
ggsave("bvardcst_gamma_s20.pdf",h)



## BVARDCST decline.length 80
h <- ggplot(data=bvardcst,aes(x=decline.length.80)) +
    geom_histogram(aes(y=..density..),fill=bvdc.col(10)[1],alpha=0.75,binwidth=max(bvardcst$decline.length.80)/100) +
        geom_density(colour=bvdc.col(10)[10],size=1.5,fill=bvdc.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Decline Length") +
                            ylab("Probability Density")
ggsave("bvardcst_decline_length_80.pdf",h)

## BVARDCST decline.length 80
h <- ggplot(data=bvardcst,aes(x=decline.length.50)) +
    geom_histogram(aes(y=..density..),fill=bvdc.col(10)[1],alpha=0.75,binwidth=max(bvardcst$decline.length.50)/100) +
        geom_density(colour=bvdc.col(10)[10],size=1.5,fill=bvdc.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Decline Length") +
                            ylab("Probability Density")
ggsave("bvardcst_decline_length_50.pdf",h)

## BVARDCST decline.length 20
h <- ggplot(data=bvardcst,aes(x=decline.length.20)) +
    geom_histogram(aes(y=..density..),fill=bvdc.col(10)[1],alpha=0.75,binwidth=max(bvardcst$decline.length.20)/100) +
        geom_density(colour=bvdc.col(10)[10],size=1.5,fill=bvdc.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Decline Length") +
                            ylab("Probability Density")
ggsave("bvardcst_decline_length_20.pdf",h)












## BCSTDVAR ntip 80
h <- ggplot(data=bcstdvar,aes(x=ntip80)) +
    geom_histogram(aes(y=..density..),fill=bcdv.col(10)[1],alpha=0.75,binwidth=max(bcstdvar$ntip80)/100) +
        geom_density(colour=bcdv.col(10)[10],size=1.5,fill=bcdv.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Number of Tips") +
                            ylab("Probability Density")
ggsave("bcstdvar_ntip_80.pdf",h)

## BCSTDVAR ntip 50
h <- ggplot(data=bcstdvar,aes(x=ntip50)) +
    geom_histogram(aes(y=..density..),fill=bcdv.col(10)[1],alpha=0.75,binwidth=max(bcstdvar$ntip50)/100) +
        geom_density(colour=bcdv.col(10)[10],size=1.5,fill=bcdv.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Number of Tips") +
                            ylab("Probability Density")
ggsave("bcstdvar_ntip_50.pdf",h)

## BCSTDVAR ntip 20
h <- ggplot(data=bcstdvar,aes(x=ntip20)) +
    geom_histogram(aes(y=..density..),fill=bcdv.col(10)[1],alpha=0.75,binwidth=max(bcstdvar$ntip20)/100) +
        geom_density(colour=bcdv.col(10)[10],size=1.5,fill=bcdv.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Number of Tips") +
                            ylab("Probability Density")
ggsave("bcstdvar_ntip_20.pdf",h)

## BCSTDVAR ntip short20
h <- ggplot(data=bcstdvar,aes(x=ntip_s20)) +
    geom_histogram(aes(y=..density..),fill=bcdv.col(10)[1],alpha=0.75,binwidth=max(bcstdvar$ntip_s20)/100) +
        geom_density(colour=bcdv.col(10)[10],size=1.5,fill=bcdv.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Number of Tips") +
                            ylab("Probability Density")
ggsave("bcstdvar_ntip_s20.pdf",h)



## BCSTDVAR gamma 80
h <- ggplot(data=bcstdvar,aes(x=gamma.80)) +
    geom_histogram(aes(y=..density..),fill=bcdv.col(10)[1],alpha=0.75) +
        geom_density(colour=bcdv.col(10)[10],size=1.5,fill=bcdv.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Corrected Gamma") +
                            ylab("Probability Density")
ggsave("bcstdvar_gamma_80.pdf",h)

## BCSTDVAR gamma 80
h <- ggplot(data=bcstdvar,aes(x=gamma.50)) +
    geom_histogram(aes(y=..density..),fill=bcdv.col(10)[1],alpha=0.75) +
        geom_density(colour=bcdv.col(10)[10],size=1.5,fill=bcdv.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Corrected Gamma") +
                            ylab("Probability Density")
ggsave("bcstdvar_gamma_50.pdf",h)

## BCSTDVAR gamma 20
h <- ggplot(data=bcstdvar,aes(x=gamma.20)) +
    geom_histogram(aes(y=..density..),fill=bcdv.col(10)[1],alpha=0.75) +
        geom_density(colour=bcdv.col(10)[10],size=1.5,fill=bcdv.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Corrected Gamma") +
                            ylab("Probability Density")
ggsave("bcstdvar_gamma_20.pdf",h)

## BCSTDVAR gamma short20
h <- ggplot(data=bcstdvar,aes(x=gamma.s20)) +
    geom_histogram(aes(y=..density..),fill=bcdv.col(10)[1],alpha=0.75) +
        geom_density(colour=bcdv.col(10)[10],size=1.5,fill=bcdv.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Corrected Gamma") +
                            ylab("Probability Density")
ggsave("bcstdvar_gamma_s20.pdf",h)




## BCSTDVAR decline.length 80
h <- ggplot(data=bcstdvar,aes(x=decline.length.80)) +
    geom_histogram(aes(y=..density..),fill=bcdv.col(10)[1],alpha=0.75,binwidth=max(bcstdvar$decline.length.80)/100) +
        geom_density(colour=bcdv.col(10)[10],size=1.5,fill=bcdv.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Corrected Decline.Length") +
                            ylab("Probability Density")
ggsave("bcstdvar_decline_length_80.pdf",h)

## BCSTDVAR decline.length 80
h <- ggplot(data=bcstdvar,aes(x=decline.length.50)) +
    geom_histogram(aes(y=..density..),fill=bcdv.col(10)[1],alpha=0.75,binwidth=max(bcstdvar$decline.length.50)/100) +
        geom_density(colour=bcdv.col(10)[10],size=1.5,fill=bcdv.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Corrected Decline.Length") +
                            ylab("Probability Density")
ggsave("bcstdvar_decline_length_50.pdf",h)

## BCSTDVAR decline.length 20
h <- ggplot(data=bcstdvar,aes(x=decline.length.20)) +
    geom_histogram(aes(y=..density..),fill=bcdv.col(10)[1],alpha=0.75,binwidth=max(bcstdvar$decline.length.20)/100) +
        geom_density(colour=bcdv.col(10)[10],size=1.5,fill=bcdv.col(10)[10],alpha=0.25) +
            theme_bw() +
#                xlim(0,10) +
#                    ylim(0,0.25) +
                        xlab("Corrected Decline.Length") +
                            ylab("Probability Density")
ggsave("bcstdvar_decline_length_20.pdf",h)



## Calculating gamma "normal"
bvdc.trees80 <- lapply(as.list(paste0("~/Dropbox/bamm_panda/bvardcst/trees/80/",1:2000,"_mol80.txt")),read.tree)
bvdc.trees50 <- lapply(as.list(paste0("~/Dropbox/bamm_panda/bvardcst/trees/50/",1:2000,"_mol50.txt")),read.tree)
bvdc.trees20 <- lapply(as.list(paste0("~/Dropbox/bamm_panda/bvardcst/trees/20/",1:2000,"_mol20.txt")),read.tree)
bvdc.trees.s20 <- lapply(as.list(paste0("~/Dropbox/bamm_panda/bvardcst/trees/short20/",1:2000,"_s20mol.txt")),read.tree)

bcdv.trees80 <- lapply(as.list(paste0("~/Dropbox/bamm_panda/bcstdvar/trees/80/",1:2000,"_80mol.txt")),read.tree)
bcdv.trees50 <- lapply(as.list(paste0("~/Dropbox/bamm_panda/bcstdvar/trees/50/",1:2000,"_50mol.txt")),read.tree)
bcdv.trees20 <- lapply(as.list(paste0("~/Dropbox/bamm_panda/bcstdvar/trees/20/",1:2000,"_20mol.txt")),read.tree)
bcdv.trees.s20 <- lapply(as.list(paste0("~/Dropbox/bamm_panda/bcstdvar/trees/short20/",1:2000,"_short20mol.txt")),read.tree)

gamma.bvdc.80 <- unlist(lapply(bvdc.trees80,gammaStat))
gamma.bvdc.50 <- unlist(lapply(bvdc.trees50,gammaStat))
gamma.bvdc.20 <- unlist(lapply(bvdc.trees20,gammaStat))
gamma.bvdc.s20 <- unlist(lapply(bvdc.trees.s20,gammaStat))

gamma.bcdv.80 <- unlist(lapply(bcdv.trees80,gammaStat))
gamma.bcdv.50 <- unlist(lapply(bcdv.trees50,gammaStat))
gamma.bcdv.20 <- unlist(lapply(bcdv.trees20,gammaStat))
gamma.bcdv.s20 <- unlist(lapply(bcdv.trees.s20,gammaStat))

bvardcst$gamma.orig.80 <- gamma.bvdc.80
bvardcst$gamma.orig.50 <- gamma.bvdc.50
bvardcst$gamma.orig.20 <- gamma.bvdc.20
bvardcst$gamma.orig.s20 <- gamma.bvdc.s20

bcstdvar$gamma.orig.80 <- gamma.bcdv.80
bcstdvar$gamma.orig.50 <- gamma.bcdv.50
bcstdvar$gamma.orig.20 <- gamma.bcdv.20
bcstdvar$gamma.orig.s20 <- gamma.bcdv.s20

write.table(bvardcst, "./bvardcst.csv", quote = FALSE, row.names = FALSE, sep = ",")
write.table(bcstdvar, "./bcstdvar.csv", quote = FALSE, row.names = FALSE, sep = ",")


#### Violin plots

fulldata <- data.frame("gamma.c80" = c(bvardcst$gamma.80,bcstdvar$gamma.80), "gamma.c50" = c(bvardcst$gamma.50,bcstdvar$gamma.50), "gamma.c20" = c(bvardcst$gamma.20,bcstdvar$gamma.20), "gamma.cs20" = c(bvardcst$gamma.s20,bcstdvar$gamma.s20), "ntip80" = c(bvardcst$ntip.80,bcstdvar$ntip80), "ntip50" = c(bvardcst$ntip.50,bcstdvar$ntip50), "ntip20" = c(bvardcst$ntip.20,bcstdvar$ntip20), "ntips20" = c(bvardcst$ntip.s20,bcstdvar$ntip_s20), "declen80" = c(bvardcst$decline.length.80, bcstdvar$decline.length.80), "declen50" = c(bvardcst$decline.length.50, bcstdvar$decline.length.50), "declen20" = c(bvardcst$decline.length.20, bcstdvar$decline.length.20))

fulldata <- data.frame("gamma.c" = c(bvardcst$gamma.80, bvardcst$gamma.50, bvardcst$gamma.20, bvardcst$gamma.s20, bcstdvar$gamma.80, bcstdvar$gamma.50, bcstdvar$gamma.20, bcstdvar$gamma.s20), "ntip" = c(bvardcst$ntip.80, bvardcst$ntip.50, bvardcst$ntip.20, bvardcst$ntip.s20, bcstdvar$ntip80, bcstdvar$ntip50, bcstdvar$ntip20, bcstdvar$ntip_s20), "declen" = c(bvardcst$decline.length.80, bvardcst$decline.length.50, bvardcst$decline.length.20, rep(NA,2000), bcstdvar$decline.length.80, bcstdvar$decline.length.50, bcstdvar$decline.length.20, rep(NA,2000)), "gamma" = c(gamma.bvdc.80, gamma.bvdc.50, gamma.bvdc.20, gamma.bvdc.s20, gamma.bcdv.80, gamma.bcdv.50, gamma.bcdv.20, gamma.bcdv.s20), "time" = factor(rep(rep(c("-80%","-50%","-20%","short20"),each=2000),2), levels = c("short20", "-20%", "-50%", "-80%")), "scenario" = c(rep("sp_var",8000),rep("ex_var",8000)))

ggplot(data=fulldata, aes(x=1, colour = scenario, fill = scenario)) +
    geom_violin(aes(y=gamma.c), scale = "width", trim = FALSE) +
    scale_colour_manual(values=c(bvdc.col(10)[10],bcdv.col(10)[10]), guide = FALSE) +
    scale_fill_manual(values=c(paste0(bvdc.col(10)[5],"90"),paste0(bcdv.col(10)[5],"90")), name = "Scenario") +
    #geom_point(aes(y=median(gamma.c), colour = scenario)) +
    theme_bw() +
    xlab("Time Slice") +
    ylab("Corrected Gamma") +
    theme(legend.position = "bottom", axis.ticks = element_blank(), axis.text.x = element_blank()) +
    facet_grid(.~time)

ggsave("gamma_corrected_per_time.pdf")

ggplot(data=fulldata, aes(x = 1, colour = scenario, fill = scenario)) +
    geom_violin(aes(y=gamma), scale = "width", trim = FALSE) +
    scale_colour_manual(values=c(bvdc.col(10)[10],bcdv.col(10)[10]), guide = FALSE) +
    scale_fill_manual(values=c(paste0(bvdc.col(10)[5],"90"),paste0(bcdv.col(10)[5],"90")), name = "Scenario") +
    theme_bw() +
    xlab("Time Slice") +
    ylab("Decline Length (MY)") +
    theme(legend.position = "bottom", axis.ticks = element_blank(), axis.text.x = element_blank()) +
    facet_grid(.~time)

ggsave("gamma_per_time.pdf")

ggplot(data=fulldata, aes(x = 1, colour = scenario, fill = scenario)) +
    geom_violin(aes(y=ntip), scale = "width", trim = FALSE) +
    scale_colour_manual(values=c(bvdc.col(10)[10],bcdv.col(10)[10]), guide = FALSE) +
    scale_fill_manual(values=c(paste0(bvdc.col(10)[5],"90"),paste0(bcdv.col(10)[5],"90")), name = "Scenario") +
    theme_bw() +
    xlab("Time Slice") +
    ylab("Number of Tips") +
    theme(legend.position = "bottom", axis.ticks = element_blank(), axis.text.x = element_blank()) +
    facet_grid(.~time)

ggsave("ntip_per_time.pdf")

ggplot(data=fulldata, aes(x = 1, colour = scenario, fill = scenario)) +
    geom_violin(aes(y=declen), scale = "width", trim = FALSE) +
    scale_colour_manual(values=c(bvdc.col(10)[10],bcdv.col(10)[10]), guide = FALSE) +
    scale_fill_manual(values=c(paste0(bvdc.col(10)[5],"90"),paste0(bcdv.col(10)[5],"90")), name = "Scenario") +
    theme_bw() +
    xlab("Time Slice") +
    ylab("Decline Length (MY)") +
    theme(legend.position = "bottom", axis.ticks = element_blank(), axis.text.x = element_blank()) +
    facet_grid(.~time)

ggsave("decline_length_per_time.pdf")





