"NewdataCom" <-
function(data1.all,data2.all,num.ori,num.class)
{
  new.data1.all <- vector("list",num.ori)
  new.data2.all <- vector("list",num.ori)  
  for ( l in 1:num.ori ){
      
      data1 <- as.matrix(data1.all[[l]])
      if (num.class == 2) {data2 <- as.matrix(data2.all[[l]])}

      temp.data <- Newdata(data1,data2,num.class)
      new.data1.all[[l]] <- temp.data$new.data1
      new.data2.all[[l]] <- temp.data$new.data2
  }
  if(num.class == 1) {new.data2.all <- NULL} 
  ##force new.data2.all to be NULl in one-sample case
  list(new.data1.all = new.data1.all,new.data2.all = new.data2.all)
}
