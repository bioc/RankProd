"RP.advance" <- function(data,cl,origin, logged=TRUE,na.rm=TRUE,gene.names=NULL,
                        plot=FALSE, rand = NULL,
                        calculateProduct=TRUE, MinNumOfValidPairs=NA,
                        RandomPairs=NA, huge=FALSE, fast=TRUE, tail.time=0.05){
# the paramenter rand is used to put the random number generator in a
# reproducible state using the rand value as seed 
    if (!is.null(rand)) { set.seed(rand)} 
    if (is.vector(data)) {
        cat("Warning: There is only one gene in the data set","\n","\n")
        data=matrix(data,nrow=1)  ##force into matrix format
    }
if(is.data.frame(data)){data <- as.matrix(data)}
check.cl <- sort(unique(cl), decreasing = F)
if (all(check.cl!= c(0,1))){
    cat("\n the class label has to contain only 0 or 1")
    stop()
}
# number of pairs equal to the number of the smaples in the smallest class
    npairs <- min(length(which(cl==0)),length(which(cl==1)))
    if (npairs==0){
        nrep<-dim(data)[2]  #in the paired case, npairs=0 but I still need to
                            #know the number of replicates for the evaluation
                            #of the pvalues
    }else{
    nrep<-npairs
    }
## Number of RandomPairs to generate, if set to NA, the odd integer
##closer to n*n is used
    if (is.na(RandomPairs)) RandomPairs <- npairs^2
    if (RandomPairs==0) RandomPairs=1
    if ((RandomPairs %% 2) == 0) RandomPairs <- RandomPairs + 1
# Replace -Inf, Inf, NaN with NA.
    data <- mzmatch.R.wrongvalues(data)
## Estimate detection limit value.
detection.limit <- min(data[data!=0],na.rm=TRUE)
#if logged=FALSE zeros are substitute with the detection.limit
if (logged==FALSE) data[data==0]<-detection.limit
## check number of samples 
    total.sam=dim(data)[2]
    total.sam1=length(cl)
    total.sam2=length(origin)
    if(total.sam!=total.sam1 | total.sam!=total.sam2)
    stop("Number of classes and/or origins should match
    the columns in the data")
# I have to separate the dataset according to the different origin
    data.pre=OriginxyCall(data,cl,origin, !calculateProduct)
#data.pre$data1 is a list of matrix, each matrix has a different origin and all
#the samples in it are from class1
#data.pre$data2 is a list of matrix, each matrix has a different origin and all
#the samples in it are from class2
    num.ori=length(unique(origin)) #counting the number of different origins
    num.gene=dim(data)[1]  #defining the number of variables
##set seed for random number generator
####AveFC
    yy<-data.pre$data2
    if (!is.null(yy) ){   ##two-sample
        num.class=2
        data1.all=data.pre$data1 ##data under condition1
        data2.all=data.pre$data2  ##data under condition2 
        fold.change=NULL
        for (l in 1:num.ori){
            data1=as.matrix(data1.all[[l]])
            data2=as.matrix(data2.all[[l]])
            data1.ave=apply(data1,1,mean) 
            data2.ave=apply(data2,1,mean) 
            if (logged) {
                fold.change=cbind(fold.change,(data1.ave-data2.ave))
            }else{
                fold.change=cbind(fold.change,(data1.ave/data2.ave))}
                rm(data1,data2,data1.ave,data2.ave)
            }
        ave.fold.change=apply(fold.change,1,mean)
        }
if ( is.null(yy)){ ##one-sample
    num.class=1
    data1.all=data.pre$data1 ##data under condition1
    data2.all=data.pre$data2  ##data under condition2 
    fold.change=NULL
    for (l in 1:num.ori){
        data1=as.matrix(data1.all[[l]])
        fold.change=cbind(fold.change,apply(data1,1,mean))
        rm(data1)
    }
    ave.fold.change=apply(fold.change,1,mean)
    }
########AveFC
#This two variables are used to save the ranks of each iteration,
#if the parameter huge is set to true these variables 
#are not saved in order to save memory
ranks1<- vector("list",RandomPairs)
ranks2<- vector("list",RandomPairs)
if(!huge){
    Orirank<-vector("list",2)
    Orirank.ori<-vector("list",2)#
    } 
badrows<-NULL
prod_out <- list ()
prod_out.Class2 <- list()
output<-vector("list", RandomPairs)
output.Class2<-vector("list", RandomPairs)
for (ori in 1:num.ori){
    data.ori <-cbind(data.pre$data1[[ori]],data.pre$data2[[ori]])
    if(!is.null(data.pre$data2)){
        cl.ori <-c(rep(0,dim(data.pre$data1[[ori]])[2]),
                    rep(1,dim(data.pre$data2[[ori]])[2]))
    }else{
        cl.ori<-c(rep(0,dim(data.pre$data1[[ori]])[2]))
    }
    npairs.ori <- min(length(which(cl.ori==0)),length(which(cl.ori==1)))
    if(npairs.ori==0) {nrep.ori <- length(cl.ori)} else {nrep.ori <- npairs.ori}
    if (is.na(MinNumOfValidPairs)) MinNumOfValidPairs <- floor(nrep.ori/2)
    xy.out.ori <- xyCall(cl.ori, RandomPairs, !calculateProduct)
    cl.ori <- xy.out.ori$cl
    x.ori <- xy.out.ori$x
    y.ori <- xy.out.ori$y
# Replace -Inf, Inf, NaN with NA.
#data.ori <- mzmatch.R.wrongvalues(data.ori)
## Estimate detection limit value.
detection.limit.ori <- min(data.ori[data.ori!=0],na.rm=TRUE)
## Prepare ratios for paired samples
    gr1.ori <- data.ori[, x.ori]
    gr2.ori <- data.ori[, y.ori]
## If in one group there was not detection in all samples,
## but there are signals detected in group B,
## values are replaced with detection.limit value (are we sure of this?)
## Also check which metabolites doesn not match MinNumOfValidPairs criteria. 
## Variables gr1NA and gr2NA will store information about which measurements
## don't contain values and those will be ignored during random pairing.
    badrows.ori <- NULL
    gr1NA.ori <- vector ("list", nrow(gr1.ori))
    gr2NA.ori <- vector ("list", nrow(gr2.ori))
    if (RandomPairs>1){
        for (coln in 1:nrow(gr1.ori)) {
            toreplacegr1.ori <- which(is.na(gr1.ori[coln,]))
            toreplacegr2.ori <- which(is.na(gr2.ori[coln,]))
            gr1NA.ori[[coln]] <- toreplacegr1.ori
            gr2NA.ori[[coln]] <- toreplacegr2.ori
            toreplacegr1.ori <- length(toreplacegr1.ori)
            toreplacegr2.ori <- length(toreplacegr2.ori)
            if (MinNumOfValidPairs > min((c(ncol(gr1.ori)-toreplacegr1.ori,
                    ncol(gr2.ori)-toreplacegr2.ori)))){
                badrows.ori <- append(badrows.ori,coln)
                gr1NA.ori[[coln]] <- NULL
                gr2NA.ori[[coln]] <- NA
            }
        if (toreplacegr1.ori==ncol(gr1.ori)){
        toreplacegr1.ori <- 1
        }else{
        toreplacegr1.ori <- 0
        }
        if (toreplacegr2.ori==ncol(gr2.ori)){
        toreplacegr2.ori <- 1
        }else{
        toreplacegr2.ori <- 0
        }
        if (toreplacegr1.ori==1 & toreplacegr2.ori==0){
        gr1.ori[coln,] <- detection.limit.ori
        }
        if (toreplacegr2.ori==1 & toreplacegr1.ori==0){
            gr2.ori[coln,] <- detection.limit.ori
        }
    }
}
    
# If option is set all NA left are replaced with the median value of the
# single feature.
    if (na.rm==TRUE){
        tmp<-NAreplace(gr1.ori,gr2.ori,y.ori)
        gr1.ori<-tmp[[1]]
        gr2.ori<-tmp[[2]]
        rm(tmp)
    }
## Create list of random sample matching, and ratios. If RandomPairs is
## set to 1, original sample order is used.
## If 'RandomPairs' is set to 1, only initial sample order is used for
## pairing (A paired test).
    ratios.list.ori <- vector("list",RandomPairs)
    ratios.list.Class2.ori <- vector("list",RandomPairs)
    for(i in 1:RandomPairs){
        ratios.list.Class2.ori[[i]] <- matrix(data = NA, ncol = npairs.ori,
                                                nrow = nrow(gr1.ori))
    }
    for(i in 1:RandomPairs){
        ratios.list.ori[[i]] <- matrix(data = NA, ncol = npairs.ori,
                                        nrow = nrow(gr1.ori))
    }
    if (is.null(y.ori)){
        ratios.list.ori[[1]] <- gr1.ori
        fold.change.ori <- apply(gr1.ori,1,mean,na.rm=TRUE)
        ratios.list.Class2.ori[[1]] <- 1/gr1.ori
        ratios.list.Class2.ori[[1]] <-
                            mzmatch.R.wrongvalues(ratios.list.Class2.ori[[1]])
    }else{
        data1.ave.ori=apply(gr1.ori, 1, mean, na.rm=TRUE)
        data2.ave.ori=apply(gr2.ori, 1, mean, na.rm=TRUE)
    if(logged==TRUE){
        ratios.list.ori[[1]] <- gr1.ori[,1:npairs.ori] - gr2.ori[,1:npairs.ori]
        ratios.list.ori[[1]] <- mzmatch.R.wrongvalues(ratios.list.ori[[1]])
        fold.change.ori=data1.ave.ori-data2.ave.ori
    }else{
        ratios.list.ori[[1]] <- gr1.ori[,1:npairs.ori]/gr2.ori[,1:npairs.ori]
        ratios.list.ori[[1]] <- mzmatch.R.wrongvalues(ratios.list.ori[[1]])
        fold.change.ori=data1.ave.ori/data2.ave.ori
    }
    if(logged==TRUE){
    ratios.list.Class2.ori[[1]] <- gr2.ori[,1:npairs.ori]-gr1.ori[,1:npairs.ori]
    ratios.list.Class2.ori[[1]] <-
                            mzmatch.R.wrongvalues(ratios.list.Class2.ori[[1]])
    }else{
    ratios.list.Class2.ori[[1]] <- gr2.ori[,1:npairs.ori]/gr1.ori[,1:npairs.ori]
    ratios.list.Class2.ori[[1]] <-
                            mzmatch.R.wrongvalues(ratios.list.Class2.ori[[1]])
        }  
    }
 
if (RandomPairs>1){
    for (coln in 1:nrow(gr1.ori)){
        val1.ori <- gr1.ori[coln,]
        val2.ori <- gr2.ori[coln,]
        ## get rid of NA's
        val1.ori <- sort(val1.ori)
        val2.ori <- sort(val2.ori)
        # calculate ratios of ranodmly paired samples
        for (rp in 2:RandomPairs) {
            val1gr.ori <- sample(val1.ori)[1:npairs.ori]
            val2gr.ori <- sample(val2.ori)[1:npairs.ori]
            if (logged==TRUE){
                ratios.list.ori[[rp]][coln,] <- val1gr.ori-val2gr.ori
                ratios.list.Class2.ori[[rp]][coln,] <- val2gr.ori-val1gr.ori
            }else{
                ratios.list.ori[[rp]][coln,] <- val1gr.ori/val2gr.ori
                ratios.list.Class2.ori[[rp]][coln,] <- val2gr.ori/val1gr.ori
            }
        }
    }
}
badrows<-unique(append(badrows,badrows.ori))
       
###########################################################################
## Ranks calculation
for (randpair in 1:RandomPairs){
    if (logged==TRUE){
        inputmatrix.ori <- ratios.list.ori[[randpair]]
    }else{
        MIN <- min(ratios.list.ori[[randpair]],na.rm=TRUE)
        if (MIN <=0) {
            cat ("Error: Input data contain negative values or values equal
                to 0.","\n")
            stop ()  
        }else{
            inputmatrix.ori <- log (ratios.list.ori[[randpair]], 10)
        }
    }
      
## replace values in badrows with NA, after "order" function rank to those rows
## will be set to NA, and ignored in sum or prod calculation
    if (length(badrows.ori)!=0){
        inputmatrix.ori[badrows.ori,] <- NA
      }
output.ori <- matrix(ncol=ncol(inputmatrix.ori),nrow=nrow(inputmatrix.ori))
inputmatrix.Class2.ori <- inputmatrix.ori*(-1)
output.Class2.ori <- matrix(ncol=ncol(inputmatrix.Class2.ori),
                            nrow=nrow(inputmatrix.Class2.ori))

for (coln in 1:ncol(inputmatrix.ori)){
# After assigning rank with "order" function, ranks for NA's are set back to NA
    NAs <- which(is.na(inputmatrix.ori[,coln]))
    if (length(NAs)>0){
        output.ori[order(inputmatrix.ori[,coln],na.last=NA),coln] <-
            c(1:(nrow(inputmatrix.ori)-length(NAs)))
        }else{
            output.ori[order(inputmatrix.ori[,coln]),coln] <-
                c(1:nrow(inputmatrix.ori))
        }
    }
for (coln in 1:ncol(inputmatrix.Class2.ori)){
    NAs <- which(is.na(inputmatrix.Class2.ori[,coln]))
    if (length(NAs)>0){
    output.Class2.ori[order(inputmatrix.Class2.ori[,coln],na.last=NA),coln]<-
    c(1:(nrow(inputmatrix.Class2.ori)-length(NAs)))
        }else{
            output.Class2.ori[order(inputmatrix.Class2.ori[,coln]),coln] <-
                                c(1:nrow(inputmatrix.Class2.ori))
        }
    }
#saving the ranks of each pairing, if the parameter huge is set to true
#these variables are not saved in order to save memory
if(ori>1){
    ranks1[[randpair]]<-cbind(output.ori,ranks1[[randpair]])
    ranks2[[randpair]]<-cbind(output.Class2.ori,ranks2[[randpair]])
}else{
    ranks1[[randpair]]<-output.ori
    ranks2[[randpair]]<-output.Class2.ori
    }
}
if (!is.null(badrows)){
    cat("Warning: There are",length(badrows),"variables with at least one
        missing value. You can consider using a 'MinNumOfValidPairs' of this
        function.","\n","\n")
    cat(" These values won't be used to compute rank product.","\n","\n")
}

########################################################################
########################################################################
#### evaluating the orirank output as it was in the old version##
if(!huge & RandomPairs>1){
    Orirank.ori<-computeOrirank(gr1.ori,gr2.ori,logged)
    Orirank[[1]]<-cbind(Orirank[[1]],Orirank.ori[[1]])
    Orirank[[2]]<-cbind(Orirank[[2]],Orirank.ori[[2]])
}

}
for( jj in 1:RandomPairs){
# performing the Rank Product or the Rank Sum according
# to the parameter calculateProduct
output<-ranks1[[jj]]
output.Class2<-ranks2[[jj]]
    if (calculateProduct==TRUE){
        ProdFunc <- function(x, DATA) exp(mean(log(DATA[x,]),na.rm=TRUE))
        products <- sapply (1:nrow(output), ProdFunc, DATA=output)
        products.Class2 <- sapply (1:nrow(output.Class2),ProdFunc,
                                   DATA=output.Class2)
    } else {
        products <- apply(output,1,mean,na.rm=TRUE)
        products.Class2 <- apply (output.Class2,1, mean,na.rm=TRUE)
    }
prod_out[[jj]] <- products
prod_out.Class2[[jj]] <- products.Class2
}
prod_out <- do.call (rbind,prod_out)
ranks_out <- apply(prod_out,2,median,na.rm=TRUE)
prod_out.Class2 <- do.call (rbind,prod_out.Class2)
ranks_out.Class2 <- apply(prod_out.Class2,2,median,na.rm=TRUE)
### evaluation of pvalues and pfps
# if we are performing the RankProduct we use the fast algorithm to evaluate
# the pfps
pvals<-matrix(NA,nrow=nrow(data), ncol=2)
pfps<-matrix(NA,nrow=nrow(data),ncol=2)
if (calculateProduct==TRUE){
    ###check max value
    MaxRP <- which(ranks_out>nrow(data))
    MaxRP2 <- which(ranks_out.Class2>nrow(data))
    ranks_out[MaxRP] <- nrow(data)
    ranks_out.Class2[MaxRP2] <- nrow(data)
    #
    ind<-which(!is.na(ranks_out))
    try(pvals[ind,1]<-rankprodbounds(ranks_out[ind]^nrep,
            (nrow(data)-length(badrows)),nrep,Delta='geometric'),silent=TRUE)
    if (is.na(pvals[ind[1],1])){
        for (i in 1:nrow(pvals)){
            if(is.na(ranks_out[i])){
                pvals[i,1]<-NA
            }else {
                pvals[i,1]<-rankprodbounds(ranks_out[i]^nrep,
                            (nrow(data)-length(badrows)),nrep,Delta='geometric')
            }
        }
    }
ind<-which(!is.na(ranks_out.Class2))
try(pvals[ind,2]<-rankprodbounds(ranks_out.Class2[ind]^nrep,
            (nrow(data)-length(badrows)),nrep,Delta='geometric'),silent=TRUE)
if (is.na(pvals[ind[1],2])){
    for (i in 1:nrow(pvals)){
        if(is.na(ranks_out.Class2[i])){
            pvals[i,2]<-NA
        }else {
            pvals[i,2]<-rankprodbounds(ranks_out.Class2[i]^nrep,
                    (nrow(data)-length(badrows)),nrep,Delta='geometric')
        }
    }
}
}else{ ### pvalues evaluation for Rank Sums
    ## using the fast function to evaluate the pvalues
    ij <- which(!is.na(ranks_out))
    ji <- which(!is.na(ranks_out.Class2))
    if(((nrep <= 30) |
        (nrep>=30 & nrep < 40 & (nrow(data)-length(badrows)<= 10^7)) |
        (nrep>=40 & nrep <= 50 & (nrow(data)-length(badrows)<= 10^6)) |
        (nrep>=50 & nrep <= 60 & (nrow(data)-length(badrows)<= 5*10^4)) |
        (nrep>=60 & nrep <= 70 & (nrow(data)-length(badrows)<= 10^4)))){
        
        pvals[,1]<- sapply(ranks_out[ij]*nrep, dice.sum.cdf, n=nrep,
                           s=nrow(data)-length(badrows))
        
        pvals[,2]<- sapply(ranks_out.Class2[ji]*nrep, dice.sum.cdf, n=nrep,
                           s=nrow(data)-length(badrows))
        
    }else if(fast == FALSE){
        
        pvals[,1]<- sapply(ranks_out[ij]*nrep, dice.sum.cdf.accurate, n=nrep,
                           s=nrow(data)-length(badrows))
        
        pvals[,2]<- sapply(ranks_out.Class2[ji]*nrep, dice.sum.cdf.accurate, n=nrep,
                           s=nrow(data)-length(badrows))
    } else{
        ord.RS1 <- order(ranks_out, na.last = TRUE)
        ord.RS2 <- order(ranks_out.Class2, na.last = TRUE)
        jhi <-1
        time <- 0
        str <- Sys.time()
        while((time < tail.time) & (jhi < (nrow(data)-length(badrows)))){
            pvals[ord.RS1[jhi],1]<-
                dice.sum.cdf.accurate(ranks_out[ord.RS1[jhi]]*nrep,
                                      nrep,(nrow(data)-length(badrows)))
            pvals[ord.RS2[jhi],2]<-
                dice.sum.cdf.accurate(ranks_out.Class2[ord.RS2[jhi]]*nrep,
                                      nrep,(nrow(data)-length(badrows)))  
            jhi<-jhi+1
            time <- as.numeric(Sys.time()-str, units="mins")
        }
        cat("\nThe highest p-values evaluated with the exact method are",
            pvals[ord.RS1[jhi-1],1],"and", pvals[ord.RS1[jhi-1],2])
        NA1 <- which(is.na(pvals[,1]))
        NA2 <- which(is.na(pvals[,2]))
        pvals[NA1,1]<-pnorm(ranks_out[NA1]*nrep,nrep*(((nrow(data)-length(badrows))+1)/2),
                               sqrt(((((nrow(data)-length(badrows))^2)-1)/12)*nrep))
        pvals[NA2,2]<-pnorm(ranks_out.Class2[NA2]*nrep,nrep*(((nrow(data)-length(badrows))+1)/2),
                               sqrt(((((nrow(data)-length(badrows))^2)-1)/12)*nrep))
    }
}
### evaluation of the pfps from the pvalues
rankorder1 <- rep(NA,nrow(pvals))
rankorder1[order(pvals[,1],na.last=NA)] <- 1:length(order(pvals[,1],na.last=NA))
rankorder2 <- rep(NA,nrow(pvals))
rankorder2[order(pvals[,2],na.last=NA)] <- 1:length(order(pvals[,2],na.last=NA))
for (i in 1:nrow(pfps)){
    pfps[i,1] <- (pvals[i,1]*(nrow(data)-length(badrows)))/rankorder1[i]
    pfps[i,2] <- (pvals[i,2]*(nrow(data)-length(badrows)))/rankorder2[i]
}  
# output RPs (RSs) and RPrank (RSrank)
out <- list()
if (calculateProduct==TRUE){
    out$RPs <- cbind(ranks_out,ranks_out.Class2)
    RPrank<-matrix(NA,nrow(data),2)
    RPrank[,1]<-rank(out$RPs[,1], na.last =TRUE)
    RPrank[badrows,1]<-NA
    RPrank[,2]<-rank(out$RPs[,2], na.last =TRUE)
    RPrank[badrows,2] <- NA
    colnames (RPrank) <- c("class1 < class2","c11lass1 > class2")
    colnames(out$RPs) <- c("class1 < class2","class1 > class2")
    out$RPrank<-RPrank
} else {
    out$RSs <- cbind(ranks_out,ranks_out.Class2)
    RSrank<-matrix(NA,nrow(data),2)
    RSrank[,1]<-rank(out$RSs[,1], na.last =TRUE)
    RSrank[badrows,1]<-NA
    RSrank[,2]<-rank(out$RSs[,2], na.last =TRUE)
    RSrank[badrows,2] <- NA
    colnames (RSrank) <- c("class1 < class2","class1 > class2")
    colnames(out$RSs) <- c("class1 < class2","class1 > class2")
    out$RSrank<-RSrank
}
# output pfps, pvals, AveFC, groups, RandomPairs, nrep
colnames (pfps) <- c("class1 < class2","class1 > class2")
colnames (pvals) <- c("class1 < class2","class1 > class2") 
out$pfp<-pfps
out$pval<-pvals
fold.change=t(t(ave.fold.change))
fold.change[badrows]<-NA
out$AveFC<-fold.change
out$groups <- cl
out$RandomPairs_ranks <- prod_out
out$nrep <- nrep
# if the parameter huge is set to TRUE the allrank1, allrank2 and
# Orirank are not saved in order to save space!
if(!huge){
    if(RandomPairs >1){
        allrank1<-ranks1
        allrank2<-ranks2
        names(allrank1)=c("class1 < class2","class1 > class 2")
        names(allrank2)=c("class1 < class2","class1 > class 2")
        out$allrank1<-allrank1
        out$allrank2<-allrank2
        out$Orirank<-Orirank
#### evaluating the orirank output as it was in the old version##
        out$Orirank<-Orirank
    }
}  
if(!is.null(gene.names)){
    out<-genenames(out,gene.names,RandomPairs)
}
if(plot){
    plotPFP(out,calculateProduct)
}
cat ("\n done  ")  
out
}
