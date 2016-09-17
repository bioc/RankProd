"NAreplace"<-function(gr1,gr2,y){
for (hh in 1:dim(gr1)[1]){
    for (gg in 1:dim(gr1)[2]){
        if (is.na(gr1[hh,gg])){
            gr1[hh,gg] <- median(c(gr1[hh,],gr2[hh,]), na.rm=TRUE)
        }
    }
    if( !is.null(y) ){      
        for (gg in 1:dim(gr2)[2]){
            if (is.na(gr2[hh,gg])){
                gr2[hh,gg] <- median(c(gr1[hh,],gr2[hh,]), na.rm=TRUE)
            }
        }
    }
    }
out<-vector("list",2)
out[[1]] <- gr1
if( !is.null(y) ) out[[2]] <- gr2
out
}