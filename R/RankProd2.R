"RankProd2" <-
function(data1.all,data2.all,num.ori,num.gene,logged,num.class,rev.sorting)
{ 
   num.rank.all=0
   rank.rep.all=t(t(1:num.gene)) ##column vector
   for (l in 1:num.ori)
    { data1=data1.all[[l]]
      data2=data2.all[[l]]

      data1=as.matrix(data1)
      if (num.class==2) {data2=as.matrix(data2)}  ##in case of one replication

      temp=RankComp(data1,data2,logged,num.class,rev.sorting)
      rank.rep.all=cbind(rank.rep.all,temp$rank)
      num.rank.all=num.rank.all+temp$num.rank
      rm(temp)
     }
 
    rank.all=rank.rep.all[,-1]

    ##compute rank product
    rank.prod.temp=rank.all^(1/num.rank.all)
    rank.prod=apply(rank.prod.temp,1,prod)
    #rank.prod=(apply(rank.all,1,prod))^(1/num.rank.all)
    rank.prod[num.rank.all==0]=NA ## genes with all value missing, assigned to NA

    list(RP=rank.prod,rank.all=rank.all)
}
