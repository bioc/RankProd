"dice.sum.cdf.accurate" <- function (p, n, s) {
    # p is the number of total points obtained
    # n is the number of dice
    # s is the number of sides on each die
    a <- trunc(n*(s+1)/2)
    if(p<=a){
        P <- 0
        for (i in 0:floor((p-n)/s)){
            P <- P + (-1)^i * chooseZ(n,i) * chooseZ(p-(s*i),n)
        }
    }else{
        P <- 0
        for (i in 0:floor((s*n-p)/s)){
            P <- P + (-1)^i * chooseZ(n,i) * chooseZ(n*(s+1)-p-s*i-1, n)
        }
    }
    if(p<=a){
        P<- P/mpfr(pow.bigz(s,n), precBits=2048)
    }else{
        P<- 1- P/mpfr(pow.bigz(s,n), precBits=2048)
    }
    return(as.numeric(P))
}