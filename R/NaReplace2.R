"NaReplace2" <-
function(X,origin){
    ori.lev <- unique(origin)
    uni.ori <- length(ori.lev)   
    for(i in 1:nrow(X)) {
        for ( c in 1:uni.ori) {       
           X.i <- X[i,origin == ori.lev[c]]
           X.i <- replace(X.i, which(is.na(X.i)), mean(X.i, na.rm = TRUE))
           X[i,origin == ori.lev[c]] <- X.i
           rm(X.i)
        }
    }
    return(X)
}
