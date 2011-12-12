"RankProd2V2" <-
function(data1.all,data2.all,num.ori,num.gene,logged,num.class,rev.sorting)
{ 
   num.rank.all=0
   rank.rep.all=NULL ## fb: this now stores a list of FUNCTION objects
   num.reps=NULL

   for (l in 1:num.ori)
    { data1=data1.all[[l]]
      data2=data2.all[[l]]

      data1=as.matrix(data1)
      if (num.class==2) {data2=as.matrix(data2)}  ##in case of one replication

      temp=RankCompV2(data1,data2,logged,num.class,rev.sorting)
      rank.rep.all=c(rank.rep.all,temp$rank) 
      num.rank.all=num.rank.all+temp$num.rank
      num.reps = cbind(num.reps, temp$num.rep)
      rm(temp)
     }
 
#    rank.all=rank.rep.all[,-1] ##TODO what about this
    rank.all <- "removed due to optimization"  ## fb: we do not have this object
    rank.prod <- rep(1,num.gene);	## fb: initialize ranks

    ##compute rank product
## fb: Use the function objects to do this without large matrices in memory
    for (part in 1:length(num.reps)) {
      num.col <- num.reps[part];
      rank.rep.getCol <- rank.rep.all[[part]];  ## fb: get the function object providing columns
      for(col in 1:num.col) {
        rank.prod.temp <- rank.rep.getCol(col)^(1/num.rank.all);  
        rank.prod <- rank.prod * rank.prod.temp;
      }
#    rank.prod.temp=rank.all^(1/num.rank.all)
#    rank.prod=apply(rank.prod.temp,1,prod)
    }
##
    #rank.prod=(apply(rank.all,1,prod))^(1/num.rank.all)
    rank.prod[num.rank.all==0]=NA ## genes with all value missing, assigned to NA

    list(RP=rank.prod,rank.all=rank.all)
}
