data <- read.csv("./results/fulldata_bvardcst.csv")

tree <- read.tree("~/Box Sync/bammpanda_timetravel/cloud/trees/80/complete/1070_comp.txt")
tree.mol <- read.tree("./trees/80/1070_mol80_9.588_0.997_0.173_1.409_104_.txt")

source("~/Dropbox/collaborations/Le Codes/branching.times.with.extinction.R")


pitoco <- data$time.80[1070]-max(branching.times.with.extinct(tree))


# Function to calculate species richness in a given point in time
rich <- function(x,parms){
    a <- parms[1]
    b <- parms[2]
    c <- parms[3]
    res <- exp((a*(-exp(-b*x))+a-(b*c*x))/b);return(res)
}



par(mfrow=c(2,1))


pdf("phylo_timetravel.pdf")
plot(tree,show.tip.label=FALSE);axisPhylo()
abline(v=max(branching.times.with.extinct(tree)),col=2, lty=2, lwd=2)
abline(v=data$time.50[1070], col=3, lty=2, lwd=2)
abline(v=data$time.20[1070], col=4, lty=2, lwd=2)
abline(v=data$time.s20[1070], col=6, lty=2, lwd=2)
dev.off()

plot(y=rich(x=seq(pitoco,data$time.80[1070],0.01)-pitoco,parms=c(data$lambda.init.sim.80[1070],data$alpha.sim.80[1070],data$mu.sim.80[1070])),x=seq(pitoco,data$time.80[1070],0.01)-pitoco,type="l",ylab="Species Richness",xlab="Time")
abline(v=max(branching.times.with.extinct(tree)),col=2, lty=2)
abline(v=data$time.50[1070], col=3, lty=2)
abline(v=data$time.20[1070], col=4, lty=2)
abline(v=data$time.s20[1070], col=6, lty=2)


rtt <- data.frame("rich"=rich(x=seq(pitoco,data$time.80[1070],0.01)-pitoco,parms=c(data$lambda.init.sim.80[1070],data$alpha.sim.80[1070],data$mu.sim.80[1070])),"time"=seq(pitoco,data$time.80[1070],0.01)-pitoco)

library(ggplot2)
library(devtools)
source("http://bioconductor.org/biocLite.R")
biocLite("ggtree")
library(ggtree)

g1 <- ggplot(data=rtt,aes(x=time,y=rich)) + geom_line() + theme_bw() + geom_vline(xintercept=max(branching.times.with.extinct(tree)),linetype="longdash",colour="red") + geom_vline(xintercept=data$time.50[1070],linetype="longdash",colour="green3") + geom_vline(xintercept=data$time.20[1070],linetype="longdash",colour="blue") + geom_vline(xintercept=data$time.s20[1070],linetype="longdash", colour="magenta") + xlim(0,4) + annotate("text",x=max(branching.times.with.extinct(tree)) + 0.2,y=250,label="80%",colour="red") + annotate("text",x=data$time.50[1070] - 0.2,y=250,label="50%",colour="green3") + annotate("text",x=data$time.20[1070] - 0.2,y=250,label="20%",colour="blue") + annotate("text",x=data$time.s20[1070] - 0.35,y=250,label="short20%",colour="magenta") + theme(axis.text.x=element_blank()) + ylab("Species richness") + xlab("Time")

ggsave("diversity_timetravel.pdf",g1)



ggsave("ntip80.pdf",ggplot(data=data,aes(x=ntip.80)) + geom_histogram(binwidth=1) + theme_bw() + xlab("Number of Species") + ylab("Count"))
ggsave("simtime.pdf",ggplot(data=data,aes(x=time.80)) + geom_histogram(binwidth=1) + theme_bw() + xlab("Tree Length (Myrs)") + ylab("Count"))
ggsave("gamma80.pdf",ggplot(data=data,aes(x=gamma.80)) + geom_histogram(binwidth=0.01) + theme_bw() + geom_vline(xintercept=0,linetype="longdash",colour="red") + xlab("Corrected Gamma Statistic") + ylab("Count"))
