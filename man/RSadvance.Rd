\name{RSadvance}
\alias{RSadvance}
\title{Advanced Rank Sum Analysis of Microarray}
\description{
    Advance rank sum method to identify differentially 
     expressed genes. It is possible to combine data from 
   different studies, e.g. data sets generated at different 
   laboratories.}

\usage{
    RSadvance(data,cl,origin,num.perm=100,logged=TRUE,
                 na.rm=FALSE,gene.names=NULL,plot=FALSE, rand=NULL)
}

\arguments{
      \item{data}{the data set that should be analyzed. Every 
      row of this data set must correspond to a gene.}

      \item{cl}{a vector containing the class labels of the 
               samples. In the two class unpaired case, the label 
               of a sample is either 0 (e.g., control group) or 1 
               (e.g., case group). For one group data, the label for 
               each sample should be 1.}

    \item{origin}{a vector containing the origin labels of the 
              sample. e.g. for 
            the data sets generated at multiple laboratories, the label
            is the same for samples within one lab and different for samples 
            from different labs. } 

  \item{num.perm}{number of permutations used in the calculation 
               of the null density. Default is 'B=100'.}
  

   \item{logged}{if "TRUE", data has bee logged, otherwise set 
                it to "FALSE"}


    \item{na.rm}{if 'FALSE' (default), the NA value will not
             be used in computing rank. If 'TRUE', the missing 
            values will be replaced by the genewise mean of
             the non-missing values. Gene will all value missing 
             will be assigned "NA"}

    \item{gene.names}{if "NULL", no gene name will be attached 
                to the estimated percentage of false prediction (pfp). }

     \item{plot}{If "TRUE", plot the estimated pfp verse the rank 
                of each gene}

      \item{rand}{if specified, the random number generator 
               will be put in a  reproducible state.}
}

\value{
      A result of identifying differentially expressed 
     genes between two classes. The identification consists of two parts,
      the identification of  up-regulated  and down-regulated genes in class 2
      compared to class 1, respectively. 

     \item{pfp}{estimated percentage of false positive predictions
             (pfp) up to  the position of each gene under two 
              identificaiton each}
     \item{pval}{estimated pvalue for each gene being up- and down-regulated}  

     \item{RSs}{Origina rank-sum (average rank) of each genes} 
    
     \item{RSrank}{rank of the rank sum of each gene in ascending order}

     \item{Orirank}{original ranks in each comparison, which is 
                 used to compute rank sum}

     \item{AveFC}{ fold change of average expression under class 1 over 
                 that under class 2, if multiple origin, than avraged 
                 across all origin. log-fold change if data is in log scaled, 
                 original fold change if data is unlogged. }
     \item{all.FC}{fold change of class 1/class 2 under each origin.
                   log-fold change if data is in log scaled}
}

\note{Percentage of false prediction (pfp), in theory, is 
      equivalent of false  discovery rate (FDR), and it is 
      possible to be large than 1.
 
      The function looks for up- and down- regulated genes in two
     seperate steps, thus two pfps are computed and used to identify 
     gene that belong to each group. 
  
     The function is able to deal with single or multiple-orgin 
      studies. It is similar to  funcion RP.advance expect a rank
      sum is computed instead of rank product. This method 
     is more sensitive to individual rank values, while rank 
     product is more robust to 
     outliers (refer RankProd vignette for details) 
}


\author{Fangxin Hong \email{fhong@salk.edu}}

\seealso{
    \code{\link{topGene}}   \code{\link{RP}}  
     \code{\link{plotRP}}  \code{\link{RPadvance}}
}



\examples{
      
      #Suppose we want to check the consistence of the data 
      #sets generated in two different 
      #labs. For example, we would look for genes that were \
      # measured to be up-regulated in 
      #class 2 at lab 1, but down-regulated in class 2 at lab 2.\
       data(arab)
      arab.cl2 <- arab.cl

      arab.cl2[arab.cl==0 &arab.origin==2] <- 1

      arab.cl2[arab.cl==1 &arab.origin==2] <- 0

      arab.cl2
  ##[1] 0 0 0 1 1 1 1 1 0 0


      #look for genes differentially expressed
      #between hypothetical class 1 and 2
      arab.sub=arab[1:500,] ##using subset for fast computation
      arab.gnames.sub=arab.gnames[1:500]
      Rsum.adv.out <- RSadvance(arab.sub,arab.cl2,arab.origin,
                          num.perm=100,
logged=TRUE,
                          gene.names=arab.gnames.sub,rand=123)

      attributes(Rsum.adv.out)
      
}

\keyword{htest}

