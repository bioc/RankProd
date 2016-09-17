"computeOrirank"<-function(gr1,gr2,logged){
    k1=dim(gr1)[2]
    k2=dim(gr2)[2]
    num.rep=k1*k2
    data.rep1=matrix(NA, nrow(gr1),num.rep)
if (logged) {
    for (k in 1:k1){
        temp=((k-1)*k2+1):(k*k2)
        data.rep1[,temp]=gr1[,k]-gr2
    }
    } else { 
        for (k in 1:k1) {
            temp=((k-1)*k2+1):(k*k2)
            data.rep1[,temp]=gr1[,k]/gr2 
            }   
    }
temp<-gr1
gr1<-gr2
gr2<-temp
rm(temp)
k1=dim(gr1)[2]
k2=dim(gr2)[2]
num.rep=k1*k2
data.rep2=matrix(NA, nrow(gr1),num.rep)
if (logged) {
    for (k in 1:k1){
        temp=((k-1)*k2+1):(k*k2)
        data.rep2[,temp]=gr1[,k]-gr2
    }
    } else { 
        for (k in 1:k1) {
            temp=((k-1)*k2+1):(k*k2)
            data.rep2[,temp]=gr1[,k]/gr2 
        }   
    }

rank.rep1=apply(data.rep1,2,rank)  ##rank of genes in the ascending order,
                                    ##NA on bottom
rank.rep2=apply(data.rep2,2,rank)  ##rank of genes in the ascending order,
                                    ## NA on bottom
Orirank<-list(rank.rep1,rank.rep2)
names(Orirank)=c("class1 < class2","class1 > class 2")
Orirank
}