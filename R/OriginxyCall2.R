"OriginxyCall2" <-
function (data,cl,origin)
{
   lev <- unique(cl)
   uni.cl <- length(lev)
   
   if (uni.cl > 2 )
     stop("There is errors in the classlabels")
    
   ori.lev <- unique(origin)
   uni.ori <- length(ori.lev)
   cat(" The data is from ",uni.ori,"different origins \n \n")
   if (min(ori.lev) != 1 | max(ori.lev) != uni.ori){
       cat("Warning: origins labels are not numbered from 1 to ",uni.ori,"\n","\n")}
       

   if(uni.cl == 1){
       
        data2 <- NULL
        data1 <- vector("list",uni.ori)
        for ( c in 1:uni.ori) {data1[[c]] = data[,origin == ori.lev[[c]]]}
   }
   if(uni.cl == 2) {
         
	data2 <- vector("list",uni.ori)
        data1 <- vector("list",uni.ori)
        for ( c in 1:uni.ori) {
            index1 <- which(origin == ori.lev[[c]] & cl == 0)           
            index2 <- which(origin == ori.lev[[c]] & cl == 1)
            if (length(index1) == 0 | length(index1) == 0)
                stop("Error: data from different origins should contain data from both classs")
            data1[[c]] <- data[,index1]
            data2[[c]] <- data[,index2]
            rm(index1,index2) }
        
    }
    list(data1 = data1,data2 = data2)
}
