"xyCall" <- function(cl, RandomPairs,sum=FALSE){
lev=unique(cl)
uni.cl=length(lev)
if (uni.cl>2 ) stop("There is something wrong with the classlabels")
if(uni.cl==1){
    cat("Rank Product analysis for paired case","\n","\n")
    if (lev!=1){
cat("warning: Expected class label is 1, cl will thus be set to 1.","\n","\n")
        cl=rep(1,length(cl)) 
    }
        x <- which(cl==1)
        y <- NULL
}
if(uni.cl==2){
    if (RandomPairs==1){
        if(sum){
            cat("Rank Sum analysis for paired case","\n","\n")
        }else{
            cat("Rank Product analysis for paired case","\n","\n")
        }
    }else{
        if(sum){
        cat("Rank Sum analysis for unpaired case","\n","\n")
        }else{
        cat("Rank Product analysis for unpaired case","\n","\n")
        }
    }
    if(min(lev)!=0 | max(lev)!=1){
        cat("Warning: Expected classlabels are 0 and 1. cl will thus be set
                to 0 and 1.","\n","\n")
        cl[which(cl==min(lev))]<-0 ##smallest one is assigned to 0
        cl[which(cl==max(lev))]<-1 ##largest one is assigned to 1
    }
    x<-which(cl==0) ##sample 1
    y<-which(cl==1) ##sample 2
    }
structure(list(x=x, y=y, cl=cl))
}
