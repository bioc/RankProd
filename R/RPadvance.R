"RPadvance" <-
function(data,cl,origin,num.perm=100,logged=TRUE,na.rm=FALSE,gene.names=NULL,plot=FALSE, rand=NULL, huge=FALSE)
{
  if (huge)
    return(RPadvanceV2(data,cl,origin,num.perm,logged,na.rm,gene.names,plot,rand))
  else
    return(RPadvanceV1(data,cl,origin,num.perm,logged,na.rm,gene.names,plot,rand));
 }
