"xyCall" <-
function(cl)
{
  lev=unique(cl)
  uni.cl=length(lev)
  #test=setdiff(lev,c(0,1))
  if (uni.cl>2 )
    stop("There is something wrong with the classlabels")
  if(uni.cl==1){
        cat("Rank Product analysis for one-class case","\n","\n")
        if (lev!=1) {
             cat("warning: Expected classlabel is 1, cl will hus be set to 1.","\n","\n")
             cl=rep(1,length(cl)) }
        x<-which(cl==1)
        y=NULL
   }
   if(uni.cl==2) {
         cat("Rank Product analysis for two-class case","\n","\n")
         if(min(lev)!=0 | max(lev)!=1){
			cat("Warning: Expected classlabels are 0 and 1. cl will thus be set to 0 and 1.","\n","\n")
			cl[which(cl==min(lev))]<-0 ##small one is assigned to 0
			cl[which(cl==max(lev))]<-1 ##big one is assigned to 0
		}
		x<-which(cl==0) ##sample 1
		y<-which(cl==1) ##sample 2
   }
   structure(list(x=x,y=y))
}
