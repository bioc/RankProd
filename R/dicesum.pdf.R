"dice.sum.pdf" <- function (p, n, s) {
# p is the number of total points obtained
# n is the number of dice
# s is the number of sides on each die
        i=0:floor((p-n)/s)
        P = sum(1/(s^n) * (-1)^i * choose(n,i) * choose(p-(s*i)-1, n-1) )
    return(as.numeric(P))
}