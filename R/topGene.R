
"topGene" <-
function(x,cutoff=NULL,method="pfp",num.gene=NULL,logged=TRUE,logbase=2,
         gene.names=NULL){  
   ##input is x: an RP object
   pfp=as.matrix(x$pfp)
   FC=as.matrix(x$AveFC)  ##data1/ data2
   pval=as.matrix(x$pval)
   
   if (is.null(x$RPs) ){  ##Rank Sum
      RP=as.matrix(x$RSs)
      rank=as.matrix(x$RSrank)
    } else {
      RP=as.matrix(x$RPs)
      rank=as.matrix(x$RPrank)
    }  

 
   if (is.null(num.gene) & is.null(cutoff)) 
       stop("No selection criteria is input, please input either
            cutoff or num.gene")

   
   ##for up-regulation under class2,data1<data2, two-sample
   ## under denominator class, one-sample 
   RP.sort.upin2=sort(RP[,1],index.return=TRUE)
   RP.sort.downin2=sort(RP[,2],index.return=TRUE)

   if (!is.null(cutoff) ) {

      if (method == "pfp") {
      cutgenes.upin2=which(pfp[RP.sort.upin2$ix,1]<cutoff)
      cutgenes.downin2=which(pfp[RP.sort.downin2$ix,2]<cutoff)
        } else {
          if (method == "pval") {
            cutgenes.upin2=which(pval[RP.sort.upin2$ix,1]<cutoff)
            cutgenes.downin2=which(pval[RP.sort.downin2$ix,2]<cutoff)
          } else {
          stop("No criterion is input to select genes, please select either
               pfp(fdr) or pval(P-value)")
          }
        }


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
         cat("Warning: gene.names should have the same length as
             the gene vector.","\n")
         cat("No gene.names are assigned","\n") 
         } else { 
         rownames(pfp)=gene.names
         #cat("gene.names are assigned","\n")
         } 
   }


   pfp[,2] <- as.numeric(prettyNum(pfp[,2],digits=4, with=6))
   pfp[,1] <- as.numeric(prettyNum(pfp[,1],digits=4, with=6))
   pval[,1] <- as.numeric(prettyNum(pval[,1],digits=4, with=6))
   pval[,2] <- as.numeric(prettyNum(pval[,2],digits=4, with=6))
   RP[,1] <- as.numeric(prettyNum(RP[,1],digits=4, with=6))
   RP[,2] <- as.numeric(prettyNum(RP[,2],digits=4, with=6))
   if (logged) {        
       FC[,1]=as.numeric(prettyNum(logbase^FC[,1],digits=4, width=6))
   } else { FC[,1]=as.numeric(prettyNum(FC[,1],digits=4, width=6)) }


   ##for up-regulation under class2,data1<data2, two-sample
   ## under denominator class, one-sample 
   
   if(length(gene.sel.upin2)>0) {
     Out.table.upin2=cbind(gene.sel.upin2,RP[gene.sel.upin2,1],
                           FC[gene.sel.upin2], pfp[gene.sel.upin2,1],
                           pval[gene.sel.upin2,1])
     rownames(Out.table.upin2)=rownames(pfp)[gene.sel.upin2]
     colnames(Out.table.upin2)=c("gene.index","RP/Rsum","FC:(class1/class2)",
                                 "pfp","P.value")  

     cat("Table1: Genes called significant under class1 < class2","\n\n")
    } else { 
     cat("No genes called significant under class1 < class2","\n\n") 
     Out.table.upin2=NULL
    }

   ##for down-regulation under class2, data1> data2, two-sample
   ## under numeratorclass, one-sample 
   
   if(length(gene.sel.downin2)>0) {
     Out.table.downin2=cbind(gene.sel.downin2,RP[gene.sel.downin2,2],
                             FC[gene.sel.downin2],pfp[gene.sel.downin2,2],
                             pval[gene.sel.downin2,2])
     rownames(Out.table.downin2)=rownames(pfp)[gene.sel.downin2]
     colnames(Out.table.downin2)=c("gene.index","RP/Rsum","FC:(class1/class2)",
                                   "pfp","P.value")  

     cat("Table2: Genes called significant under class1 > class2","\n\n")  
    } else { 
      cat("No genes called significant under class1 > class2","\n\n") 
      Out.table.downin2=NULL
    }


   list(Table1=Out.table.upin2,Table2=Out.table.downin2)
    
}


