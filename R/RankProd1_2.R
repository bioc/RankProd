"RankProd1v2" <-
function(data1,data2,logged,num.class,rev.sorting)
{

  ##version 2: For each Fold-change, get rank, then rankprod, only keep the rankprod updated from each Fold-change
  ##no matrix stored for memory efficiency

  ##compute all possible pairwise comparison for two sample
  ##deal with NA values (assume non-TDE value)
  ##get rank product 

  ## fb: decide on an operator for calculating fold changes
  fcoperator = if (logged) get("-") else get("/")
   

  num.gene <- dim(data1)[1]


  ##two-smaple
  if (num.class == 2) { 

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
     
     #number of rank for each gene
     num.rank <- rep(0, num.gene) 
     for ( k in 1:k1)
     {  for (j in 1:k2)
        { temp=fcoperator(data1.wk[,k], data2.wk[,j])
          num.rank <- num.rank + !is.na(temp)
        }
     }

      
     rank.prod <- rep(1,num.gene)
     for ( k in 1:k1)
       {temp=((k-1)*k2+1):(k*k2)
        temp2=fcoperator(data1.wk[,k], data2.wk)

        ##rank of genes in the ascending order, NA on bottom
        rank.temp=apply(temp2,2,rank)
        rank.temp[is.na(temp2)]=1 ## 1 won't affect the product

        rank.prod.temp=apply(rank.temp^(1/num.rank),1,prod) ##take first, also applicable to general case
        rank.prod=rank.prod *rank.prod.temp

        rm (temp,temp2,rank.temp,rank.prod.temp)
         }  

     rm(data1.wk,data2.wk,k1,k2,num.rep)

    
   }

    ##one-sample
   if ( num.class==1) { 
      # data.rep=data1 ##itself, no need
      
      if (rev.sorting) {  #rank from largest to smallest
          num.rep=dim(data1)[2]
          rank.temp=matrix(NA,num.gene,num.rep)
          for ( r in 1:num.rep) {rank.temp[,r]=rank(data1[,r],na.last = FALSE) } ##NA on top
          rank.rep=(num.gene+1)-rank.temp  ##reverse rank from largest to smallest, NA on bottom
          rm(num.rep,rank.temp)
      } else { ####rank of genes in the ascending order
          rank.rep=apply(data1,2,rank)  #NA on bottom
      }

     rank.rep[is.na(data1)]=1 ## 1 won't affect the product
     ## not memory hungry, keep it the same
     num.rank=apply(is.na(data1)==FALSE,1,sum)
     rank.prod.temp=rank.rep^(1/num.rank) ##take first, also applicable to general case
     rank.prod=apply(rank.prod.temp,1,prod)

    }  

    

    ##Finalized the rank product
          
   rank.prod[num.rank==0]=NA ## genes with all value missing, assigned to NA

   gc()

   list(RP=rank.prod)
}
