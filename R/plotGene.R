"plotGene" <-
function(gene.to.plot,x,gene.names=NULL, data,cl,origin,logged=TRUE,logbase=2)
{
   ##plot expression profile of each gene user input
   num.gene=dim(as.matrix(x$pfp))[1]

   ##names for all genes
   if (!is.null(gene.names)) {
         if (num.gene!=length(gene.names)|| length(gene.names)!= dim(data)[1] )
         {   stop("Number of genes are not consistent from inputs ") }
         
   } else {
     cat("No gene.names input","\n")
     cat("Gene names from data is used","\n") 
     gene.names=rownames(data)
   }   

   ## process data with multiple origins 
   data.pre=OriginxyCall2(data,cl,origin)
   y=data.pre$data2 ##Null for one-sample case
   num.ori=length(unique(origin))
   
   ##match tge targeting gene
   index.target=pmatch(gene.to.plot,gene.names,nomatch=NA)
   if (is.null(index.target))
       stop("The input gene is not found in the input gene names")

   ##input is x: an RP object
   pfp.all=as.matrix(x$pfp)
   pfp.all[pfp.all>1]=1
   pfp=round(pfp.all[index.target,],3)
   AveFC=round(as.matrix(x$AveFC)[index.target],3)  ##data1/ data2
   pval=round(as.matrix(x$pval)[index.target,],3)
   if ( num.ori>1 ) 
       { FC.all=round(as.matrix(x$all.FC)[index.target,],3)
   } else { FC.all=NULL}
   

   
   ##plot parameters
   pfp.string1=paste("pfp=",pfp[1],sep="")
   pfp.string2=paste("pfp=",pfp[2],sep="")
   p.string1=paste("p.val=",pval[1],sep="")
   p.string2=paste("p.val=",pval[2],sep="")
   up.string=paste(pfp.string1,p.string1,sep=" ")
   down.string=paste(pfp.string2,p.string2,sep=" ")
   
       
      
   ##plot 
   par(mar=c(4,4,6,2))
   par(mfrow=c(1,num.ori))
   for ( l in 1:num.ori)
   {
      data1=data.pre$data1[[l]][index.target,]
      data2=data.pre$data2[[l]][index.target,]
      num.sam1=length(data1)
      num.sam2=length(data2) 

      if (logged) { 
          data1=logbase^data1
          data2=logbase^data2
          FC.all2=round(logbase^FC.all,3)
          AveFC2=round(logbase^AveFC,2)  }   
     

      if ( num.sam2!=0)
      { num.sam=max(num.sam1,num.sam2,na.rm=TRUE)        
        data.l=matrix(NA,2,num.sam)
        data.l[1,1:num.sam1]=data1
        data.l[2,1:num.sam2]=data2
        rownames(data.l)=c("Class1","Class2")
        
        max.y=max(data.l,na.rm=TRUE)*1.2
        barplot(t(data.l),width = 1,space =c(0,1.2), beside = TRUE,col=c(rep("rosybrown",num.sam1),
             rep("lightgoldenrod",num.sam2)),ylim=c(0,max.y),ylab="expression",xlab=paste("data set ",l,sep=""))     
      } else {
         data.l=data1
         max.y=max(data.l,na.rm=TRUE)*1.2
         barplot(data.l,width = 1,space =0, beside = TRUE,col="burlywood2",
              ylim=c(0,max.y),ylab="expression",xlab=paste("data set ",l,sep=""))     
      }
      
        if (num.ori>1 ) 
         {text((num.sam1+3),max.y*0.9,labels=paste("F.C.=",FC.all2[l],sep=""))}
      
      if( l==1)
      {
        mtext(paste("Name:",gene.names[index.target],sep=""),line=+4,adj=0,side=3,col=1,font=3,cex=1)
        mtext(paste("Ave.FC(class1/class2)=",AveFC2,sep=""),line=+3,adj=0,side=3,col=1,font=1,cex=1)
        mtext(paste("Up-regulation (class 2):",up.string,sep=" "),line=+2,adj=0,side=3,col=1,font=1,cex=1)
        mtext(paste("down-regulation (class 2):",down.string,sep=" "),line=+1,adj=0,side=3,col=1,font=1,cex=1)
      }

      rm(data1,data2,num.sam1,num.sam2,data.l,max.y)
   }
      
  
   list(gene.name=gene.to.plot,AveFC=AveFC2, FC.all=FC.all2,pfp=pfp,pval=pval)

}
  
   
   
