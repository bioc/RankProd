"plotRP" <- function(x, cutoff = NULL){
##input is x
pfp=as.matrix(x$pfp)
if (is.null(x$RPs) ){  ##Rank Sum
    RP=as.matrix(x$RSs)
    rank=as.matrix(x$RSrank)
    } else {
        RP=as.matrix(x$RPs)
        rank=as.matrix(x$RPrank)
    }  
ind1<-which(!is.na(RP[,1]))
ind2<-which(!is.na(RP[,2]))
ind3<-append(ind1,ind2)
ind3<-unique(ind3)
##for up-regulation under class2,data1<data2, two-sample
## under denominator class, one-sample 
RP.sort.upin2=sort(RP[ind1,1],index.return=TRUE)
RP.sort.downin2=sort(RP[ind2,2],index.return=TRUE)
pfp1<-pfp[ind1,1]
pfp2<-pfp[ind2,2]
rank1<-rank[ind1,1]
rank2<-rank[ind2,2]
par(mfrow=c(2,1))
    if (!is.null(cutoff) ) {
        cutgenes.upin2=which(pfp1[RP.sort.upin2$ix]<cutoff)
        cutgenes.downin2=which(pfp2[RP.sort.downin2$ix]<cutoff)
        if (length(cutgenes.upin2) > 0) {
            numTop=max(cutgenes.upin2)  
            gene.sel2=RP.sort.upin2$ix[1:numTop]  
            plot(rank1[-gene.sel2],pfp1[-gene.sel2],type="p",
                xlim=c(1,length(pfp1)),ylim=range(pfp[ind3,]),
                xlab="number of identified genes",ylab="estimated PFP")
            title("Identification of Up-regulated genes under class 2")
            abline(h=cutoff,col="red")
            abline(v=numTop,col="red")
            par(new=TRUE)
            plot(rank1[gene.sel2],pfp1[gene.sel2],type="p",
                xlim=c(1,length(pfp1)),ylim=range(pfp[ind3,]),
                xlab="",ylab="",col="red")
            rm(numTop)
        }else {
            cat("No genes found using the input cutoff class 1 < class 2 \n")
            plot(rank1,pfp1,type="p",
                xlab="number of identified genes",ylab="estimated PFP")
            title("Identification of Up-regulated genes under class 2")
        }
        if (length(cutgenes.downin2) > 0) {
            numTop=max(cutgenes.downin2)  
            gene.sel1=RP.sort.downin2$ix[1:numTop]  
            plot(rank2[-gene.sel1],pfp2[-gene.sel1],type="p",
                xlim=c(1,length(pfp2)),ylim=range(pfp[ind3,]),
                xlab="number of identified genes",ylab="estimated PFP")
            title("Identification of down-regulated genes under class 2")
            abline(h=cutoff,col="red")
            abline(v=numTop,col="red")
            par(new=TRUE)
            plot(rank2[gene.sel1],pfp2[gene.sel1],type="p",
                xlim=c(1,length(pfp2)),ylim=range(pfp[ind3,]),
                xlab="",ylab="",col="red")
            rm(numTop)
            } else {
             cat("No genes found using the input cutoff: class 1 > class 2 \n")
             plot(rank2,pfp2,type="p",xlab="number of identified genes",
                ylab="estimated PFP")
             title("Identification of down-regulated genes under class 2")
            }
    }
if(is.null(cutoff)){
    plot(rank1,pfp1,type="p",xlim=c(1,length(pfp1)),ylim=range(pfp[ind3,]),
        xlab="number of identified genes",ylab="estimated PFP")
    title("Identification of Up-regulated genes under class 2")
    plot(rank2,pfp2,type="p",xlim=c(1,length(pfp2)),ylim=range(pfp[ind3,]),
        xlab="number of identified genes",ylab="estimated PFP")
    title("Identification of down-regulated genes under class 2")
    }
}
