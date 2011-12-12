"RPV2" <- function (data, cl, num.perm = 100, logged = TRUE, na.rm = FALSE,
    gene.names = NULL, plot = FALSE, rand = NULL)
{
    cat("Using memory-efficient version of RP\n");

    if (is.vector(data)) {
        cat("Warning: There is only one gene in the data set", "\n", "\n")
        data = matrix(data, nrow = 1)
    }

    ## fb: decide on an operator for calculating fold changes
    fcoperator = if (logged) get("-") else get("/")

    total.sam <- length(cl)
    total.sam2 <- dim(data)[2]
    if (total.sam != total.sam2)
        stop("Number of classes should match the columns in the data")

    xy.out <- xyCall(cl)
    x <- xy.out$x
    y <- xy.out$y

    num.gene <- dim(data)[1]

    ##deal with NA value
    data <- as.matrix(data)
    mode(data) <- "numeric"
    NA.genes <- NULL
    if (any(is.na(data))) { # checks if there are NAs
        NA.genes <- unique(ceiling(which(is.na(t(data)))/ncol(data)))
        cat("Warning: There are", length(NA.genes), "genes with at least one missing value.","\n", "\n")
        if (na.rm)
            data[NA.genes, ] <- NaReplace(data[NA.genes, ])
        if (!na.rm)
            cat(" This value is not used to compute rank product.",  "\n", "\n")
    }


    ##set seed for random number generator
    if (!is.null(rand)) {
        set.seed(rand)
    }

    ##make the data set ready
    if (!is.null(y)) { ##two-sample  
        num.class <- 2

        data1 <- as.matrix(data[, x])
        data2 <- as.matrix(data[, y])

        data1.ave = apply(data1, 1, mean)
        data2.ave = apply(data2, 1, mean)

        ## fb: use the fold change operator
        fold.change = fcoperator(data1.ave,data2.ave)

     } else { ##one-sample
        num.class <- 1

        data1 <- as.matrix(data[, x]) ##data under condition1
        data2 <- as.matrix(data[, y]) ##NULL  

        fold.change = apply(data1, 1, mean)
    }

    ##original rank product and rank of each rank product
    #study up-regulated genes under condition2
    RP.ori.upin2 <- RankProd1v2(data1, data2, logged, num.class, rev.sorting = FALSE)
     rank.ori.upin2 <- rank(RP.ori.upin2$RP)

    ##study down-regulated genes under condition2
    RP.ori.downin2 <- RankProd1v2(data1, data2, logged, num.class, rev.sorting = TRUE)
     rank.ori.downin2 <- rank(RP.ori.downin2$RP)

    cat("Starting", num.perm, "permutations...", "\n")
   
    ## fb: Ranks are now computed iteratively, i.e. updated with each new row of permuted data. Initialization is done here.
    temp2.up <- rank(RP.ori.upin2$RP)  ##rank.ori.upin2
    temp2.down <- rank(RP.ori.downin2$RP)  ##rank.ori.downin2

    for (p in 1:num.perm) {

        ## fb: Call the garbage collector more often
        gc()
        
        temp.data <- Newdata(data1, data2, num.class)
        new.data1 <- temp.data$new.data1
        new.data2 <- temp.data$new.data2

        ## fb: This is the new column wise iterative process: Ranks are updated right away, no permuted columns need be stored
        ##calculate RP for the permuted data
        RP.perm.upin2.oneCol <- RankProd1v2(new.data1, new.data2, logged, num.class, rev.sorting = FALSE)$RP
        RP.perm.downin2.oneCol <- RankProd1v2(new.data1, new.data2, logged, num.class, rev.sorting = TRUE)$RP
  
        
        temp2.up <- updateRank( temp2.up, RP.ori.upin2$RP, RP.perm.upin2.oneCol )
        temp2.down <- updateRank( temp2.down, RP.ori.downin2$RP, RP.perm.downin2.oneCol )

        rm(temp.data,new.data1 ,new.data2,RP.perm.upin2.oneCol,RP.perm.downin2.oneCol )

    }

    ##address significance level
    cat("Computing pfp ..", "\n")

    ##for up-regulation under class2,data1<data2, two-sample
    ## under denominator class, one-sample  
    
    temp2<-temp2.up ##the rank of original RP in the increasing order
    
    order.temp <- match(temp2, sort(temp2))
    count.perm <- (sort(temp2) - c(1:num.gene))[order.temp]
    exp.count <- count.perm/num.perm
    pval.upin2 <- count.perm/(num.perm * num.gene)
    pfp.upin2 <- exp.count/rank.ori.upin2

    if (plot) {
        par(mfrow = c(2, 1))
        plot(rank.ori.upin2, pfp.upin2, xlab = "sorted gene rank of the original rank product",
            ylab = "estimated PFP")
        title("Identification of Up-regulated genes under class 2")
    }
    rm( temp2, order.temp, count.perm, exp.count)


    ##for down-regulation under class2, data1> data2, two-sample
    ## under numeratorclass, one-sample     
   
    temp2<-temp2.down
    

    order.temp <- match(temp2, sort(temp2))
    count.perm <- (sort(temp2) - c(1:num.gene))[order.temp]
    exp.count <- count.perm/num.perm
    pval.downin2 <- count.perm/(num.perm * num.gene)
    pfp.downin2 <- exp.count/rank.ori.downin2

    if (plot) {
        plot(rank.ori.downin2, pfp.downin2, xlab = "sorted gene rank of the original rank product",
            ylab = "estimated PFP")
        title("Identification of down-regulated genes under class 2")
    }
    rm( temp2, order.temp, count.perm, exp.count)



    ##output the estimated pfp and ranks of all genes
    cat("Outputing the results ..", "\n")

    pfp = data.frame(pfp.upin2, pfp.downin2)
    colnames(pfp) = c("class1 < class2", "class1 > class 2")

    pval = data.frame(pval.upin2, pval.downin2)
    colnames(pval) = c("class1 < class2", "class1 > class 2")

    RPs = data.frame(RP.ori.upin2$RP, RP.ori.downin2$RP)
    colnames(RPs) = c("class1 < class2", "class1 > class 2")
    RPrank = data.frame(rank.ori.upin2, rank.ori.downin2)
    colnames(RPrank) = c("class1 < class2", "class1 > class 2")

    Orirank = list(RP.ori.upin2$rank.all, RP.ori.downin2$rank.all)
    names(Orirank) = c("class1 < class2", "class1 > class 2")
    fold.change = t(t(fold.change))
    colnames(fold.change) = "log/unlog(class1/class2)"

    if (!is.null(gene.names)) {
        rownames(pfp) = gene.names
        rownames(pval) = gene.names
        rownames(RPs) = gene.names
        rownames(RPrank) = gene.names
        rownames(fold.change) = gene.names
    }
    list(pfp = pfp, pval = pval, RPs = RPs, RPrank = RPrank,
        Orirank = "Orirank was not computed (huge=TRUE in call to RP)", AveFC = fold.change)
}



