"RP" <-
function(data,cl,num.perm = 100,logged = TRUE,na.rm = FALSE,gene.names = NULL,plot = FALSE, rand = NULL, huge=FALSE)
{
  if (huge)
    return(RPV2(data,cl,num.perm,logged,na.rm,gene.names,plot,rand))
  else
    return(RPV1(data,cl,num.perm,logged,na.rm,gene.names,plot,rand));
 }
