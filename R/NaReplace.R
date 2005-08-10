"NaReplace" <-
function(X){
    if (is.vector(X)) { X=matrix(X,nrow=1)}     

    for(i in 1:nrow(X))
        X[i,  ] <- replace(X[i,  ], which(is.na(X[i,  ])), mean(X[i,  ], na.rm = TRUE))
    return(X)
}
