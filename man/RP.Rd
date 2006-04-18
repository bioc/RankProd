\name{RP}
\alias{RP}
\title{Rank Product Analysis of Microarray}
\description{
   Perform rank product method to identify differentially 
   expressed genes.
   It is possible to do either a one-class or two-class analysis. 
}
\usage{
    RP(data,cl,num.perm=100,logged=TRUE,
       na.rm=FALSE,gene.names=NULL,plot=FALSE, rand=NULL)
}

\arguments{
    \item{data}{the data set that should be analyzed. Every 
           row of this data set must correspond to a gene.}
    \item{cl}{a vector containing the class labels of the samples.
            In the two class unpaired case, the label of a 
            sample is either 0 (e.g., control group) or 1 
           (e.g., case group). For one class  data, the label for 
           each sample should be 1.}
     \item{num.perm}{number of permutations used in the 
             calculation of the null density. Default is 'num.perm=100'.}
     
    \item{logged}{if "TRUE", data has bee logged, otherwise set it 
               to "FALSE"}

     \item{na.rm}{if 'FALSE' (default), the NA value will not
               be used in computing rank. If 'TRUE', the missing 
              values will be replaced by the gene-wise mean of
             the non-missing values. Gene with all values missing 
                will be assigned "NA"}
    \item{gene.names}{if "NULL", no gene name will be assigned 
                     to the estimated percentage of 
                 false positive predictions (pfp).} 
     \item{plot}{If "TRUE", plot the estimated pfp verse the 
                  rank of each gene.}
     \item{rand}{if specified, the random number generator will 
           be put in a reproducible state using the rand value as seed.}
          
}


\value{A result of identifying differentially expressed genes 
      between two classes. The identification consists of two parts,
      the identification of  up-regulated  and down-regulated genes in 
      class 2 compared to class 1, respectively.
        

     \item{pfp}{estimated percentage of false positive predictions
             (pfp) up to  the position of each gene under two 
              identificaiton each} 
     \item{pval}{estimated pvalue for each gene being up- and down-regulated} 
                 
     \item{RPs}{Original rank-product of each genes for two 
                dentificaiton each }
     \item{RPrank}{rank of the rank product of each genes}
     \item{Orirank}{original rank in each comparison, which 
                    is used to construct rank product}
     \item{AveFC}{ fold change of average expression under class 1 over 
                 that under class 2. log-fold change if data is in log 
                 scaled, original fold change if data is unlogged. }
                 
}

\note{Percentage of false prediction (pfp), in theory, is 
      equivalent of false 
     discovery rate (FDR), and it is possible to be large than 1.
  
     The function looks for up- and down- regulated genes in two
     seperate steps, thus two pfps and pvalues are computed and used to identify 
     gene that belong to each group.   

     This function is suitable to deal with data from a 
     single origin, e.g. single  experiment. If the data has 
      different origin, e.g. generated at different 
     laboratories, please refer RP.advance.
}

\author{Fangxin Hong \email{fhong@salk.edu}}

\seealso{
    \code{\link{topGene}}   \code{\link{RPadvance}}  
     \code{\link{plotRP}}
}


\references{
      Breitling, R., Armengaud, P., Amtmann, A., and Herzyk, 
       P.(2004) Rank Products:A simple, yet powerful, new method to 
        detect differentially regulated genes in
      replicated microarray experiments, \emph{FEBS Letter}, 57383-92
}

\examples{     
      # Load the data of Golub et al. (1999). data(golub) 
      # contains a 3051x38 gene expression
      # matrix called golub, a vector of length called golub.cl 
      # that consists of the 38 class labels,
      # and a matrix called golub.gnames whose third column 
      # contains the gene names.
      data(golub)

 
      #use a subset of data as example, apply the rank 
      #product method
      subset <- c(1:4,28:30)
      #Setting rand=123, to make the results reproducible,

      RP.out <- RP(golub[,subset],golub.cl[subset],rand=123) 
      
      # class 2: label =1, class 1: label = 0
      #pfp for identifying genes that are up-regulated in class 2 
      #pfp for identifying genes that are down-regulated in class 2 
      head(RP.out$pfp)
  

}

\keyword{htest}

