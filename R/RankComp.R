"RankComp" <-
function(data1,data2,logged,num.class,rev.sorting)
{
  ##compute all possible pairwise comparison for two sample
  ##deal with NA values (assume non-TDE value)
  ##get rand product 

  num.gene=dim(data1)[1]

  if (num.class==2) { ##two-smaple

     if (rev.sorting) { ##switch data1 and data2)        
        data1.wk <- data2
        data2.wk <-data1        
      }  else {
        data1.wk <- data1
        data2.wk <- data2
     } 

     k1=dim(data1.wk)[2]
     k2=dim(data2.wk)[2]
     num.rep=k1*k2
     data.rep=matrix(NA, num.gene,num.rep)

     if (logged) {
       for ( k in 1:k1)
       {temp=((k-1)*k2+1):(k*k2)
        data.rep[,temp]=data1.wk[,k]-data2.wk}
     } else { 
       for ( k in 1:k1) 
        { temp=((k-1)*k2+1):(k*k2)
          data.rep[,temp]=data1.wk[,k]/data2.wk }
     }   
     
     rank.rep=apply(data.rep,2,rank)  ##rank of genes in the ascending order, NA on bottom
   }

   if ( num.class==1) {  ##one-sample
      data.rep=data1 ##itself
      
      if (rev.sorting) {  #rank from largest to smallest
          num.rep=dim(data1)[2]
          rank.temp=matrix(NA,num.gene,num.rep)
          for ( r in 1:num.rep) {rank.temp[,r]=rank(data1[,r],na.last = FALSE) } ##NA on top
          rank.rep=(num.gene+1)-rank.temp  ##reverse rank from largest to smallest, NA on bottom
      } else { ####rank of genes in the ascending order
          rank.rep=apply(data1,2,rank)  #NA on bottom
      }
    }  

    ##compute the product
   rank.rep[is.na(data.rep)]=1 ## 1 won't affect the product
   num.rank=apply(is.na(data.rep)==FALSE,1,sum) ## no of rank each gene has

   list(rank=rank.rep,num.rank=num.rank)
}
