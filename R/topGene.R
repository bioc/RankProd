"topGene" <-
function(x,cutoff=NULL,num.gene=NULL,logged=TRUE,logbase=2,gene.names=NULL)
{  
   ##input is x: an RP object
   pfp=as.matrix(x$pfp)
   FC=as.matrix(x$AveFC)  ##data1/ data2
   
   if (is.null(x$RPs) ){  ##Rank Sum
      RP=as.matrix(x$RSs)
      rank=as.matrix(x$RSrank)
    } else {
      RP=as.matrix(x$RPs)
      rank=as.matrix(x$RPrank)
    }  

 
   if (is.null(num.gene) & is.null(cutoff)) 
       stop("No selection criteria is input, please input either cutoff or num.gene")

   
   ##for up-regulation under class2,data1<data2, two-sample
   ## under denominator class, one-sample 
   RP.sort.upin2=sort(RP[,1],index.return=TRUE)
   RP.sort.downin2=sort(RP[,2],index.return=TRUE)

   if (!is.null(cutoff) ) {
      cutgenes.upin2=which(pfp[RP.sort.upin2$ix,1]<cutoff)
      cutgenes.downin2=which(pfp[RP.sort.downin2$ix,2]<cutoff)
      
      if (length(cutgenes.upin2)>0) {
              numTop=max(cutgenes.upin2)  
              gene.sel.upin2=RP.sort.upin2$ix[1:numTop]  
              rm(numTop)
        }else {
		gene.sel.upin2=c()
       }
         
      if (length(cutgenes.downin2)>0) {
              numTop=max(cutgenes.downin2)  
              gene.sel.downin2=RP.sort.downin2$ix[1:numTop]  
              rm(numTop)
        }else {
		gene.sel.downin2=c()
       }

   }
       
   if (is.null(cutoff) & !is.null(num.gene)) {     
     if (num.gene>0) {
		gene.sel.upin2=RP.sort.upin2$ix[1:num.gene] 
                gene.sel.downin2=RP.sort.downin2$ix[1:num.gene] 

	} else {
		gene.sel.upin2=c()
                gene.sel.downin2=c()
	}
   }

   if (!is.null(gene.names)) {
         if (dim(pfp)[1]!=length(gene.names) ){
         cat("Warning: gene.names should have the same length as the gene vector.","\n")
         cat("No gene.names are assigned","\n") 
         } else { 
         rownames(pfp)=gene.names
         #cat("gene.names are assigned","\n")
         } 
   }


   pfp=round(pfp,4)
   RP=round(RP,4)
   if (logged) {        
       FC=round(logbase^FC,4)
   } else { FC=round(FC,4) }


   ##for up-regulation under class2,data1<data2, two-sample
   ## under denominator class, one-sample 
   
   if(length(gene.sel.upin2)>0) {
     Out.table.upin2=cbind(gene.sel.upin2,RP[gene.sel.upin2,1],FC[gene.sel.upin2],
                   pfp[gene.sel.upin2,1])
     rownames(Out.table.upin2)=rownames(pfp)[gene.sel.upin2]
     colnames(Out.table.upin2)=c("gene.index","RP/Rsum","FC:(class1/class2)","pfp")  

     cat("Table1: Genes called significant under class1 < class2","\n\n")
    } else { 
     cat("No genes called significant under class1 < class2","\n\n") 
     Out.table.upin2=NULL
    }

   ##for down-regulation under class2, data1> data2, two-sample
   ## under numeratorclass, one-sample 
   
   if(length(gene.sel.downin2)>0) {
     Out.table.downin2=cbind(gene.sel.downin2,RP[gene.sel.downin2,2],FC[gene.sel.downin2],
                   pfp[gene.sel.downin2,2])
     rownames(Out.table.downin2)=rownames(pfp)[gene.sel.downin2]
     colnames(Out.table.downin2)=c("gene.index","RP/Rsum","FC:(class1/class2)","pfp")  

     cat("Table2: Genes called significant under class1 > class2","\n\n")  
    } else { 
      cat("No genes called significant under class1 > class2","\n\n") 
      Out.table.downin2=NULL
    }


   list(Table1=Out.table.upin2,Table2=Out.table.downin2)
    
}


