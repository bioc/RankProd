"Newdata" <-
function(data1,data2,num.class)
{
  
  k1 <- dim(data1)[2]
  num.gene <- dim(data1)[1]
  new.data1 <- matrix(NA,num.gene,k1)
  
  for (k_1 in 1:k1)
  {
    temp <- sample(1:num.gene,num.gene,replace = FALSE, prob = NULL)
    new.data1[,k_1] <- data1[temp,k_1];
    rm(temp)
   }
   rm(k_1,k1)

  if (num.class == 2) {
     k2 <- dim(data2)[2]
     new.data2 <- matrix(NA,num.gene,k2) 
     for (k_2 in 1:k2) {   
         temp <- sample(1:num.gene,num.gene,replace = FALSE, prob = NULL)
         new.data2[,k_2] <- data2[temp,k_2];
     }
      rm(k_2,k2)
   }

   if (num.class == 1) { new.data2=NULL}
   list(new.data1 = new.data1,new.data2 = new.data2)
   
 }
