"genenames"<-function(out,gene.names,RandomPairs){
    if(!is.null(out$RPs)) rownames(out$RPs)=gene.names
    if(!is.null(out$RPrank)) rownames(out$RPrank)=gene.names
    if(!is.null(out$RSs)) rownames(out$RSs)=gene.names  
    if(!is.null(out$RSrank)) rownames(out$RSrank)=gene.names 
    if(!is.null(out$Orirank[[1]])) rownames(out$Orirank[[1]])=gene.names
    if(!is.null(out$Orirank[[2]])) rownames(out$Orirank[[2]])=gene.names
    if(!is.null(out$pfp)) rownames(out$pfp)=gene.names  
    if(!is.null(out$pval)) rownames(out$pval)=gene.names
    if(!is.null(out$AveFC)) rownames(out$AveFC)=gene.names
    if(!is.null(out$allrank1)){
    for (ih in 1:RandomPairs){
        rownames(out$allrank1[[ih]])=gene.names
        rownames(out$allrank2[[ih]])=gene.names
        }
    }
out
}



