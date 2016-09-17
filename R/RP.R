"RP" <- function(data, cl, num.perm = 100, logged = TRUE, na.rm = TRUE,
                gene.names = NULL,plot = FALSE, rand = NULL, huge=FALSE){
return(RankProducts(data,cl,logged,na.rm,gene.names,plot,rand,
        calculateProduct=TRUE, MinNumOfValidPairs=NA, RandomPairs=NA, huge))
}
