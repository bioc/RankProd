"plotPFP"<-function(out,calculateProduct){
if(calculateProduct){
    ind1<-which(!is.na(out$RPs[,1]))
    ind2<-which(!is.na(out$RPs[,2]))
    ind3<-append(ind1,ind2)
    ind3<-unique(ind3)
    RP.sort.upin2=sort(out$RPs[ind1,1],index.return=TRUE)
    RP.sort.downin2=sort(out$RPs[ind2,2],index.return=TRUE)  
    RPupin2=rank(out$RPs[ind1,1],na.last=TRUE)
    RPdownin2=rank(out$RPs[ind2,2],na.last=TRUE)  
    par(mfrow=c(2,1))
    plot(RPupin2,out$pfp[ind1,1],
        xlab="sorted gene rank of the original rank product",
        ylab="estimated PFP",
        xlim=c(1,length(out$pfp[ind1,1])),ylim=range(out$pfp[ind3,]))
    title("Identification of Up-regulated genes under class 2")  
    plot(RPdownin2,out$pfp[ind2,2],
        xlab="sorted gene rank of the original rank product",
        ylab="estimated PFP",
        xlim=c(1,length(out$pfp[ind2,1])),ylim=range(out$pfp[ind3,]))
    title("Identification of down-regulated genes under class 2")
}else {
    ind1<-which(!is.na(out$RSs[,1]))
    ind2<-which(!is.na(out$RSs[,2]))
    ind3<-append(ind1,ind2)
    ind3<-unique(ind3)
    RS.sort.upin2=sort(out$RSs[ind1,1],index.return=TRUE)
    RS.sort.downin2=sort(out$RSs[ind2,2],index.return=TRUE)  
    RSupin2=rank(out$RSs[ind1,1],na.last=TRUE)
    RSdownin2=rank(out$RSs[ind2,2],na.last=TRUE)  
    par(mfrow=c(2,1))
    plot(RSupin2,out$pfp[ind1,1],
        xlab="sorted gene rank of the original rank product",
        ylab="estimated PFP",
        xlim=c(1,length(out$pfp[ind1,1])),ylim=range(out$pfp[ind3,]))
    title("Identification of Up-regulated genes under class 2")  
    plot(RSdownin2,out$pfp[ind2,2],
        xlab="sorted gene rank of the original rank product",
        ylab="estimated PFP",
        xlim=c(1,length(out$pfp[ind2,1])),ylim=range(out$pfp[ind3,]))
    title("Identification of down-regulated genes under class 2")
    }
}