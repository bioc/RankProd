"RSadvance" <- function(data,cl,origin,num.perm=100,logged=TRUE,na.rm=TRUE,
                    gene.names=NULL,plot=FALSE,rand=NULL, huge=FALSE,fast=TRUE,
                    tail.time=0.05){
return(RP.advance(data,cl,origin,logged,na.rm,gene.names,plot,rand,
    calculateProduct=FALSE,MinNumOfValidPairs=NA, RandomPairs=NA, huge, fast,
    tail.time))
}
 