"plotRP" <-
function(x, cutoff = NULL)
{
   ##input is x
   pfp=as.matrix(x$pfp)
   
   if (is.null(x$RPs) ){  ##Rank Sum
      RP=as.matrix(x$RSs)
      rank=as.matrix(x$RSrank)
    } else {
      RP=as.matrix(x$RPs)
      rank=as.matrix(x$RPrank)
    }  

   ##for up-regulation under class2,data1<data2, two-sample
   ## under denominator class, one-sample 
   RP.sort.upin2=sort(RP[,1],index.return=TRUE)
   RP.sort.downin2=sort(RP[,2],index.return=TRUE)

   par(mfrow=c(2,1))
   if (!is.null(cutoff) ) {
      cutgenes.upin2=which(pfp[RP.sort.upin2$ix,1]<cutoff)
      cutgenes.downin2=which(pfp[RP.sort.downin2$ix,2]<cutoff)

      if (length(cutgenes.upin2) > 0) {
              
              numTop=max(cutgenes.upin2)  
              gene.sel2=RP.sort.upin2$ix[1:numTop]  
              
              plot(rank[-gene.sel2,1],pfp[-gene.sel2,1],type="p",xlim=c(1,length(pfp[,1])),ylim=range(pfp),
                   xlab="number of identified genes",ylab="estimated PFP")
              title("Identification of Up-regulated genes under class 2")
              abline(h=cutoff,col="red")
              abline(v=numTop,col="red")
              par(new=TRUE)
              plot(rank[gene.sel2,1],pfp[gene.sel2,1],type="p",xlim=c(1,length(pfp[,1])),ylim=range(pfp),xlab="",ylab="",col="red")

              rm(numTop)
            
       } else {
             cat("No genes found using the input cutoff class 1 < class 2 \n")
             plot(rank[,1],pfp[,1],type="p",xlab="number of identified genes",ylab="estimated PFP")
             title("Identification of Up-regulated genes under class 2")
       }


       if (length(cutgenes.downin2) > 0) {
              
              numTop=max(cutgenes.downin2)  
              gene.sel1=RP.sort.downin2$ix[1:numTop]  
              
              plot(rank[-gene.sel1,2],pfp[-gene.sel1,2],type="p",xlim=c(1,length(pfp[,2])),ylim=range(pfp),
                   xlab="number of identified genes",ylab="estimated PFP")
              title("Identification of down-regulated genes under class 2")
              abline(h=cutoff,col="red")
              abline(v=numTop,col="red")
              par(new=TRUE)
              plot(rank[gene.sel1,2],pfp[gene.sel1,2],type="p",xlim=c(1,length(pfp[,2])),ylim=range(pfp),xlab="",ylab="",col="red")

              rm(numTop)
            
       } else {
             cat("No genes found using the input cutoff: class 1 > class 2 \n")
             plot(rank[,2],pfp[,2],type="p",xlab="number of identified genes",ylab="estimated PFP")
             title("Identification of down-regulated genes under class 2")
             
       }

    }

    if(is.null(cutoff) ){
       plot(rank[,1],pfp[,1],type="p",xlim=c(1,length(pfp[,1])),ylim=range(pfp),xlab="number of identified genes",ylab="estimated PFP")
       title("Identification of Up-regulated genes under class 2")
       plot(rank[,2],pfp[,2],type="p",xlim=c(1,length(pfp[,2])),ylim=range(pfp),xlab="number of identified genes",ylab="estimated PFP")
       title("Identification of down-regulated genes under class 2")
             
    }
 
   
}
