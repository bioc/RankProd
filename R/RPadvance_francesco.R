"RPadvance" <- function(data,cl,origin,num.perm=100,logged=TRUE,na.rm=TRUE,
                        gene.names=NULL,plot=FALSE, rand=NULL, huge=FALSE){
return(RP.advance(data,cl,origin,logged,na.rm,gene.names,plot,rand,
            calculateProduct=TRUE,MinNumOfValidPairs=NA, RandomPairs=NA, huge))
}
 