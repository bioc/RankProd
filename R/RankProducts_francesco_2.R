"RankProducts" <- function (data, cl, logged=TRUE, na.rm=TRUE, gene.names=NULL,
                            plot=FALSE, rand = NULL, calculateProduct=TRUE,
                            MinNumOfValidPairs=NA, RandomPairs=NA,
                            huge=FALSE,fast=TRUE, tail.time=0.05){
# the paramenter rand is used to put the random number generator in a
# reproducible state using the rand value as seed 
if (!is.null(rand)) set.seed(rand)
if (is.vector(data)) {
    cat("Warning: There is only one gene in the data set","\n","\n")
    data=matrix(data,nrow=1)  ##force into matrix format
    }
if(is.data.frame(data)){data <- as.matrix(data)}
# number of pairs equal to the number of the smaples in the smallest class
npairs <- min(length(which(cl==0)),length(which(cl==1)))
if (npairs==0){
    nrep<-dim(data)[2]  #in the paired case, npairs=0 but I still need to know
                        #the number of replicates for the evaluation of
                        #the pvalues
}else{
    nrep<-npairs
}
## Number of RandomPairs to generate, if set to NA, the odd integer closer
## to n*n is used
if (is.na(RandomPairs)) RandomPairs <- npairs^2
if (RandomPairs==0) RandomPairs=1
if ((RandomPairs %% 2) == 0) RandomPairs <- RandomPairs + 1
# MinNumOfValidPairs is a parameter that indicates the minimum number of NAs
#accepted per each variable/feature/gene
# If it is set to NA the half of the number of replicates is used
if (is.na(MinNumOfValidPairs)) MinNumOfValidPairs <- nrep/2
#This two variables are used to save the ranks of each iteration,
#if the parameter huge is set to true these variables 
#are not saved in order to save memory
if(!huge){
    ranks1<- vector("list",RandomPairs)
    ranks2<- vector("list",RandomPairs)
}
## checking if the cl vector and the data matrix are congruent
total.sam <- length(cl)
total.sam2 <- dim(data)[2]
if (total.sam != total.sam2){ 
    stop("Number of classes should match the columns in the data")
}  
## this fuction gives the index of the samples of the classes
## if y is null that mean that there is only one class (paired case)
xy.out <- xyCall(cl, RandomPairs, !calculateProduct)
cl <- xy.out$cl
x <- xy.out$x
y <- xy.out$y
# Replace -Inf, Inf, NaN with NA.
data <- mzmatch.R.wrongvalues(data)
## Estimate detection limit value.
detection.limit <- min(data[data!=0],na.rm=TRUE)
#if logged=FALSE zeros are substitute with the detection.limit
if (logged==FALSE) data[data==0]<-detection.limit
## Prepare ratios for paired samples
gr1 <- data[, x]
gr2 <- data[, y]
## If in one group there was not detection in all samples, but there are signals
## detected in group B,
## values are replaced with detection.limit value 
## Also check which metabolites doesn not match MinNumOfValidPairs criteria. 
## Variables gr1NA and gr2NA will store information about which measurements
## don't contain values and those will be ignored during random pairing.
badrows <- NULL
gr1NA <- vector ("list", nrow(gr1))
gr2NA <- vector ("list", nrow(gr2))
if (RandomPairs>1){
    for (coln in 1:nrow(gr1)){
        toreplacegr1 <- which(is.na(gr1[coln,]))
        toreplacegr2 <- which(is.na(gr2[coln,]))
        gr1NA[[coln]] <- toreplacegr1
        gr2NA[[coln]] <- toreplacegr2
        toreplacegr1 <- length(toreplacegr1)
        toreplacegr2 <- length(toreplacegr2)
        if (MinNumOfValidPairs >
            min((c(ncol(gr1)-toreplacegr1,ncol(gr2)-toreplacegr2)))){
            badrows <- append(badrows,coln)
            gr1NA[[coln]] <- NULL
            gr2NA[[coln]] <- NA
        }
        if (toreplacegr1==ncol(gr1)) toreplacegr1 <- 1 else toreplacegr1 <- 0
        if (toreplacegr2==ncol(gr2)) toreplacegr2 <- 1 else toreplacegr2 <- 0
        if (toreplacegr1==1 & toreplacegr2==0){
            gr1[coln,] <- detection.limit
        }
        if (toreplacegr2==1 & toreplacegr1==0){
            gr2[coln,] <- detection.limit
        }
    }
}
if (!is.null(badrows)){
    cat("Warning: There are",length(badrows),"variables with at least one
        missing value. You can consider using a 'MinNumOfValidPairs'
        of this function.","\n","\n")
    cat(" These values won't be used to compute rank product.","\n","\n")
}
  
# If option is set all NA left are replaced with the median value
# of the single feature.
if (na.rm==TRUE){
    tmp<-NAreplace(gr1,gr2,y)
    gr1<-tmp[[1]]
    gr2<-tmp[[2]]
    rm(tmp)
}
## Create list of random sample matching, and ratios. If RandomPairs is set
## to 1, original sample order is used. If 'RandomPairs' is set to 1, only
##initial sample order is used for pairing (A paired test).
ratios.list <- vector("list",RandomPairs)
ratios.list.Class2 <- vector("list",RandomPairs)
for(i in 1:RandomPairs){
    ratios.list.Class2[[i]] <-
        matrix(data = NA, ncol = npairs, nrow = nrow(gr1))
    }
for(i in 1:RandomPairs){
    ratios.list[[i]] <- matrix(data = NA, ncol = npairs, nrow = nrow(gr1))
    }
# If there is only one class label present (y is null), we assume that
# input is already ratios. Proceed with paired test.
# Otherwise calculate class1/class2
if (is.null(y)){
    ratios.list[[1]] <- gr1
    fold.change <- apply(gr1,1,mean,na.rm=TRUE)
    ratios.list.Class2[[1]] <- 1/gr1
    ratios.list.Class2[[1]] <- mzmatch.R.wrongvalues(ratios.list.Class2[[1]])
} else {
    data1.ave=apply(gr1, 1, mean, na.rm=TRUE)
    data2.ave=apply(gr2, 1, mean, na.rm=TRUE)
    if(logged==TRUE){
        ratios.list[[1]] <- gr1[,1:npairs] - gr2[,1:npairs]
        ratios.list[[1]] <- mzmatch.R.wrongvalues(ratios.list[[1]])
        fold.change=data1.ave-data2.ave
    }else{
        ratios.list[[1]] <- gr1[,1:npairs]/gr2[,1:npairs]
        ratios.list[[1]] <- mzmatch.R.wrongvalues(ratios.list[[1]])
        fold.change=data1.ave/data2.ave
    }
    if(logged==TRUE){
        ratios.list.Class2[[1]] <- gr2[,1:npairs]-gr1[,1:npairs]
      ratios.list.Class2[[1]] <- mzmatch.R.wrongvalues(ratios.list.Class2[[1]])
    }else{
        ratios.list.Class2[[1]] <- gr2[,1:npairs]/gr1[,1:npairs]
      ratios.list.Class2[[1]] <- mzmatch.R.wrongvalues(ratios.list.Class2[[1]])
    }  
}
if (RandomPairs>1){
    for (coln in 1:nrow(gr1)){
        val1 <- gr1[coln,]
        val2 <- gr2[coln,]
## get rid of NA's
        val1 <- sort(val1)
        val2 <- sort(val2)
# calculate ratios of ranodmly paired samples
        for (rp in 2:RandomPairs) {
            val1gr <- sample(val1)[1:npairs]
            val2gr <- sample(val2)[1:npairs]
            if (logged==TRUE){
                ratios.list[[rp]][coln,] <- val1gr-val2gr
            } else{
                ratios.list[[rp]][coln,] <- val1gr/val2gr
            }
            if (logged==TRUE){
                ratios.list.Class2[[rp]][coln,] <- val2gr-val1gr
          } else {
              ratios.list.Class2[[rp]][coln,] <- val2gr/val1gr
          }
        }
    }
}
## Rank products calculation
prod_out <- list ()
prod_out.Class2 <- list()
for (randpair in 1:RandomPairs){
    if (logged==TRUE){
        inputmatrix <- ratios.list[[randpair]]
    } else {
        MIN <- min(ratios.list[[randpair]],na.rm=TRUE)
        if (MIN <=0) {
            cat ("Error: Input data contain negative values or
                values equal to 0.","\n")
        stop ()	
        } else {
            inputmatrix <- log (ratios.list[[randpair]], 10)
        }
    }
    
## replace values in badrows with NA, after "order" function rank to those
## rows will be set to NA, and ignored in sum or prod calculation
if (length(badrows)!=0){
    inputmatrix[badrows,] <- NA
}
output <- matrix(ncol=ncol(inputmatrix),nrow=nrow(inputmatrix))
inputmatrix.Class2 <- inputmatrix*(-1)
output.Class2 <- matrix(ncol=ncol(inputmatrix.Class2),
                        nrow=nrow(inputmatrix.Class2))
for (coln in 1:ncol(inputmatrix)){
# After assigning rank with "order" function, ranks for NA's are set back to NA
    NAs <- which(is.na(inputmatrix[,coln]))
    if (length(NAs)>0){
        output[order(inputmatrix[,coln],na.last=NA),coln] <-
                    c(1:(nrow(inputmatrix)-length(NAs)))
    } else {
        output[order(inputmatrix[,coln]),coln] <- c(1:nrow(inputmatrix))
    }
}
for (coln in 1:ncol(inputmatrix.Class2)){
    NAs <- which(is.na(inputmatrix.Class2[,coln]))
        if (length(NAs)>0){
        output.Class2[order(inputmatrix.Class2[,coln],na.last=NA),coln] <-
                            c(1:(nrow(inputmatrix.Class2)-length(NAs)))
        }else{
            output.Class2[order(inputmatrix.Class2[,coln]),coln] <-
                                    c(1:nrow(inputmatrix.Class2))
        }
}
#saving the ranks of each pairing, if the parameter huge is set to true
#these variables are not saved in order to save memory
if (!huge){
    ranks1[[randpair]]<-output
    ranks2[[randpair]]<-output.Class2
}
#performing the Rank Product or the Rank Sum according to the
#parameter calculateProduct
if (calculateProduct==TRUE){
    ProdFunc <- function(x, DATA) exp(mean(log(DATA[x,]),na.rm=TRUE))
    products <- sapply (1:nrow(output),ProdFunc, DATA=output)
    products.Class2 <- sapply (1:nrow(output.Class2),ProdFunc,
                                DATA=output.Class2)
    } else {
        products <- apply(output,1,mean,na.rm=TRUE)
        products.Class2 <- apply (output.Class2,1, mean,na.rm=TRUE)
    }
prod_out[[randpair]] <- products
prod_out.Class2[[randpair]] <- products.Class2
}
prod_out <- do.call (rbind,prod_out)
ranks_out <- apply(prod_out,2,median,na.rm=TRUE)
prod_out.Class2 <- do.call (rbind,prod_out.Class2)
ranks_out.Class2 <- apply(prod_out.Class2,2,median,na.rm=TRUE)
### evaluation of pvalues and pfps
### if we are performing the RankProduct we use the
### fast algorithm to evaluate the pfps


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
    }else{ 
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
                pvals[ord.RS1[jhi-1],1],"and", pvals[ord.RS2[jhi-1],2]," ")
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
for (i in 1:nrow(data)){
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
    colnames (RPrank) <- c("class1 < class2","class1 > class2")
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
fold.change=t(t(fold.change))
fold.change[badrows]<-NA
out$AveFC<-fold.change
colnames(fold.change)="log/unlog(class1/class2)"
out$groups <- cl
out$RandomPairs_ranks <- prod_out
out$nrep <- nrep
# if the parameter huge is set to TRUE the allrank1, allrank2 and Orirank
# are not saved in order to save space
if(!huge){
    if(RandomPairs >1){
        allrank1<-ranks1
        allrank2<-ranks2
        names(allrank1)=c("class1 < class2","class1 > class 2")
        names(allrank2)=c("class1 < class2","class1 > class 2")
        out$allrank1<-allrank1
        out$allrank2<-allrank2
#### evaluating the orirank output as it was in the old version##
        out$Orirank<-computeOrirank(gr1,gr2,logged)
################################################
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
