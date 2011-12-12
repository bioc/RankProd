"RankCompV2" <-
function(data1,data2,logged,num.class,rev.sorting)
{
  ##compute all possible pairwise comparison for two sample
  ##deal with NA values (assume non-TDE value)
  ##get rand product 

## fb: decide on an operator for calculating fold changes
    fcoperator = if (logged) get("-") else get("/");
##    

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

## fb: new functions to create a vector of data on demand instead of holding the full matrix
        data.rep.getCol <- function(col) {
          col <- col - 1 # start at index zero
          i2 <- col %% k2 +1;
          i1 <- col %/% k2 +1;
          fcoperator(data1.wk[,i1], data2.wk[,i2])
        }
        rank.rep.getCol <- function(col) {
          dta <- data.rep.getCol(col)
          rnk <- rank(dta);
          rnk[is.na(dta)] <- 1;
          rnk
        }
##

## fb: compute the ranks iteratively instead of using the full matrix data.rep
        num.rank <- rep(0, num.gene);
        for (col in 1:num.rep) {
          num.rank <- num.rank + !is.na(data.rep.getCol(col))
        }
        num.rank <- as.integer(num.rank) 
##

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

## fb: The one class case is less memory hungry, so only the access to the data is made compatible to the two-class case
        data.rep.getCol <- function(col) { data1[,col] }
        rank.rep.getCol <- function(col) { rank.rep[,col] }
        num.rank <- apply(is.na(data1) == FALSE, 1, sum)
##      
    }  

    ##compute the product
## fb: rank.rep is replaced by rank.rep.getCol, num.rank has already been computed in the two cases above
#   rank.rep[is.na(data.rep)]=1 ## 1 won't affect the product
#   num.rank=apply(is.na(data.rep)==FALSE,1,sum) ## no of rank each gene has
##
   list(rank=rank.rep.getCol,num.rank=num.rank, num.rep=num.rep)  ## fb: now returning the FUNCTION for the rank columns as well as number of reps
}
