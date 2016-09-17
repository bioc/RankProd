"dice.sum.cdf" <- function (p, n, s) {
    # p is the number of total points obtained
    # n is the number of dice
    # s is the number of sides on each die
    a <- trunc(n*(s+1)/2)
    if (p <= a){
        i=0:floor((p-n)/s)
        P = sum(1/(s^n) * (-1)^i * choose(n,i) * choose(p-(s*i), n) )
        return(as.numeric(P))
    }else{
        i=0:floor((s*n-p)/s)
        P = 1-sum(1/(s^n) * (-1)^i * choose(n,i) * choose(n*(s+1)-p-s*i-1, n) )
        return(as.numeric(P))
    }
    
}
