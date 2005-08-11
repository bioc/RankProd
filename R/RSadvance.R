"RSadvance" <-
function(data,cl,origin,num.perm=100,logged=TRUE,na.rm=FALSE,gene.names=NULL,plot=FALSE, rand=NULL)
{
  ##different from RP, as it use variance as penalty

  total.sam=dim(data)[2]
  total.sam1=length(cl)
  total.sam2=length(origin)
   
  if ( total.sam!=total.sam1 | total.sam!=total.sam2)
     stop("Number of classes and/or origins should match the columns in the data")

  data.pre=OriginxyCall(data,cl,origin,sum=TRUE)
  y=data.pre$data2 ##Null for one-sample case
 
  num.ori=length(unique(origin))
  num.gene=dim(data)[1]
  
  ##deal with NA value
  data=as.matrix(data)
  mode(data)="numeric"

  NA.genes<-NULL
  if(any(is.na(data))){   # checks if there are NAs
     NA.genes<-unique(ceiling(which(is.na(t(data)))/ncol(data)))  ##row of missing   
     cat("Warning: There are",length(NA.genes),"genes with at least one missing value.","\n","\n")
     if(na.rm) 
        data[NA.genes,]<-NaReplace2(data[NA.genes,],origin) # replace missing values with the gene mean        
     if(!na.rm)
         cat(" This value is not used to compute rank product.","\n","\n")       
  }

  
  ##set seed for random number generator
  if (!is.null(rand)) { set.seed(rand)}

  ##main program, two-sample or one-sam case
  if (!is.null(y) ){   ##two-sample
      num.class=2
     
      data1.all=data.pre$data1 ##data under condition1
      data2.all=data.pre$data2  ##data under condition2 

      fold.change=matrix(0,num.gene,1)
      for (l in 1:num.ori)
      { data1=as.matrix(data1.all[[l]])
        data2=as.matrix(data2.all[[l]])
 
        data1.ave=apply(data1,1,mean) 
        data2.ave=apply(data2,1,mean) 
        if (logged) { fold.change=fold.change+(data1.ave-data2.ave)
        }else {fold.change=fold.change+(data1.ave/data2.ave)}
        rm(data1,data2,data1.ave,data2.ave)
      }
      ave.fold.change=fold.change/num.ori
      
    
   }

   
   if ( is.null(y)){ ##one-sample
       num.class=1
         
       data1.all=data.pre$data1 ##data under condition1
       data2.all=data.pre$data2  ##data under condition2 

       fold.change=matrix(0,num.gene,1)
       for (l in 1:num.ori)
       { data1=as.matrix(data1.all[[l]])
         fold.change=fold.change+apply(data1,1,mean) 
         rm(data1)
       }
       ave.fold.change=fold.change/num.ori
   }      
 
    ##original rank product and rank of each rank product
    ##study up-regulated genes under condition2 ,data1 < data2
    RS.ori.out.upin2=RankSum(data1.all,data2.all,num.ori,num.gene,logged,num.class,rev.sorting=FALSE)
    RS.ori.upin2=RS.ori.out.upin2$RS
    rank.ori.upin2=rank(RS.ori.upin2)
    ##study down-regulated genes under condition2,data1 < data2
    RS.ori.out.downin2=RankSum(data1.all,data2.all,num.ori,num.gene,logged,num.class,rev.sorting=TRUE)
    RS.ori.downin2=RS.ori.out.downin2$RS
    rank.ori.downin2=rank(RS.ori.downin2)
  
    

    ##permutation
    RS.perm.upin2=matrix(NA,num.gene,num.perm)
    RS.perm.downin2=matrix(NA,num.gene,num.perm)

    cat("Starting ",num.perm,"permutations...","\n")
    for (p in 1:num.perm) {
 
       new.data.temp=NewdataCom(data1.all,data2.all,num.ori,num.class)
       new.data1.all=new.data.temp$new.data1.all
       new.data2.all=new.data.temp$new.data2.all

       temp1=RankSum(new.data1.all,new.data2.all,num.ori,num.gene,logged,num.class,rev.sorting=FALSE)
       RS.perm.upin2[,p]=temp1$RS
       rm(temp1)
       temp2=RankSum(new.data1.all,new.data2.all,num.ori,num.gene,logged,num.class,rev.sorting=TRUE)
       RS.perm.downin2[,p]=temp2$RS
       rm(temp2)
       
    }
    
    ##address significance level
   cat("Computing pfp...","\n")

   ##for up-regulation under class2,data1<data2, two-sample
   ## under denominator class, one-sample 
   temp1=as.vector(cbind(RS.ori.upin2,RS.perm.upin2))
   temp2=rank(temp1)[1:num.gene]  ##the rank of original RP in the increasing order
   order.temp=match(temp2,sort(temp2))
   count.perm=(sort(temp2)-c(1:num.gene))[order.temp]
   exp.count=count.perm/num.perm

   
   pfp.upin2= exp.count/rank.ori.upin2
   if (plot) 
   { par(mfrow=c(2,1)) 
     plot(rank.ori.upin2,pfp.upin2,xlab="number of identified genes",ylab="estimated PFP")
     title("Identification of Up-regulated genes under class 2")
    }
   rm(temp1,temp2,order.temp,count.perm,exp.count)

   ##for down-regulation under class2, data1> data2, two-sample
   ## under numeratorclass, one-sample 
   temp1 <- as.vector(cbind(RS.ori.downin2,RS.perm.downin2))
   temp2 <- rank(temp1)[1:num.gene]  ##the rank of original RP in the increasing order
   order.temp <- match(temp2,sort(temp2))
   count.perm <- (sort(temp2)-c(1:num.gene))[order.temp]
   exp.count <- count.perm/num.perm
 
   pfp.downin2 <- exp.count/rank.ori.downin2
   if (plot) 
   { plot(rank.ori.downin2,pfp.downin2,xlab="sorted gene rank of the original rank product",ylab="estimated PFP")
     title("Identification of down-regulated genes under class 2")
    }

   rm(temp1,temp2,order.temp,count.perm,exp.count)

   ##output the estimated pfp and ranks of all genes
   pfp=data.frame(pfp.upin2,pfp.downin2)
   colnames(pfp)=c("class1 < class2","class1 > class 2")

   RSs=data.frame(RS.ori.upin2,RS.ori.downin2)
   colnames(RSs)=c("class1 < class2","class1 > class 2")

   RSrank=data.frame(rank.ori.upin2,rank.ori.downin2)
   colnames(RSrank)=c("class1 < class2","class1 > class 2")

   Orirank=list(RS.ori.out.upin2$rank.all,RS.ori.out.downin2$rank.all)
   names(Orirank)=c("class1 < class2","class1 > class 2")

   ave.fold.change=t(t(ave.fold.change))
   colnames(ave.fold.change)="log/unlog(class1/class2)"

   if (!is.null(gene.names)) {
           rownames(pfp)=gene.names
           rownames(RSs)=gene.names
           rownames(RSrank)=gene.names
           rownames(ave.fold.change)=gene.names
   } 


    
   list(pfp=pfp,RSs=RSs,RSrank=RSrank,Orirank=Orirank,AveFC=ave.fold.change)
       
    
 }
