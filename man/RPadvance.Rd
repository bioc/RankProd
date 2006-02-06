\name{RPadvance}
\alias{RPadvance}
\title{Advanced Rank Product Analysis of Microarray}
\description{Advance rank product method to identify 
      differentially expressed genes. It is possible to
    combine data from different studies, e.g. data sets 
    generated at different laboratories.}

\Usage{
    RPadvance(data,cl,origin,num.perm=100,logged=TRUE,
              na.rm=FALSE,gene.names=NULL,plot=FALSE, 
               rand=NULL)
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

\value{A result of identifying differentially expressed 
        genes between two classes. The identification consists of two parts,
      the identification of  up-regulated  and down-regulated genes in class 2
      compared to class 1, respectively. 

     \item{pfp}{estimated percentage of false positive predictions
             (pfp) up to  the position of each gene under two 
              identificaiton each}
     \item{pval}{estimated pvalue for each gene being up- and down-regulated} 
     
 
     \item{RPs}{Original rank-product of each genes for two i
                dentificaiton each } 
    
     \item{RPrank}{rank of the rank products of each gene in 
                 ascending order}

     \item{Orirank}{original ranks in each comparison, which 
                     is used to compute rank product}
     \item{AveFC}{ fold change of average expression under class 1 over 
                 that under class 2, if multiple origin, than avraged 
                 across all origin. log-fold change if data is in log scaled, 
                 original fold change if data is unlogged. }
           
}

\note{Percentage of false prediction (pfp), in theory, is 
      equivalent of false discovery rate (FDR), and it is 
      possible to be large than 1.

      The function looks for up- and down- regulated genes in two
     seperate steps, thus two pfps are computed and used to identify 
     gene that belong to each group.   

     The function is able to replace function RP in the 
     same library. it is a more  general version, as it is
        able to handle data from differnt origins.    
}

\author{Fangxin Hong \email{fhong@salk.edu}}

\seealso{
    \code{\link{topGene}}   \code{\link{RP}}  
      \code{\link{plotRP}}  \code{\link{RSadvance}}
}


\reference{
      Breitling, R., Armengaud, P., Amtmann, A., and Herzyk, 
      P.(2004) Rank Products: A simple, yet powerful, new method 
      to detect differentially regulated genes in
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

      ##For data with single origin
      subset <- c(1:4,28:30)
      origin <- rep(1,7)
      #identify genes 
      RP.out <- RPadvance(golub[,subset],golub.cl[subset],
                           origin,plot=FALSE,rand=123)
      
      #For data from multiple origins
      
      #Load the data arab in the package, which contains 
      # the expression of 22,081 genes
      # of control and treatment group from the experiments 
      #indenpently conducted at two 
      #laboratories.
      data(arab)
      arab.origin #1 1 1 1 1 1 2 2 2 2
      arab.cl #0 0 0 1 1 1 0 0 1 1
      RP.adv.out <- RPadvance(arab,arab.cl,arab.origin,
                    num.perm=100,gene.names=arab.gnames,logged=TRUE,rand=123)

      attributes(RP.adv.out)
      head(RP.adv.out$pfp)
      head(RP.adv.out$RPs)
      head(RP.adv.out$AveFC)


}


\keyword{htest}

