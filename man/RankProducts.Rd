\name{RankProducts}
\alias{RankProducts}
\title{Rank Product/Rank Sum Analysis}
\description{The function performs the Rank Product (or Rank Sum) method to
identify differentially expressed genes. It is possible to do either a
one-class or two-class analysis.}
\usage{
RankProducts(data, cl, logged = TRUE, na.rm = TRUE, gene.names = NULL,
plot = FALSE, rand = NULL, calculateProduct = TRUE, MinNumOfValidPairs = NA,
RandomPairs = NA, huge = FALSE, fast = TRUE, tail.time = 0.05)
}
\arguments{
\item{data}{the data set that should be analyzed. Every row of this dataset must
correspond to a gene}
\item{cl}{a vector containing the class labels of the samples.
In the two class unpaired case, the label of a sample
is either 0 (e.g., control group) or 1 (e.g., case group).
For one class data, the label for each sample should be 1}
\item{logged}{if "TRUE" data have been previously log transformed. Otherwise it 
should be set as "FALSE"}
\item{na.rm}{if "FALSE", the NA value will not be used in
computing rank. If "TRUE" (default), the missing values will be replaced by
the genewise median of the non-missing values.
Gene with a number of missing values greater than "MinNumOfValidPairs" are still
not considered in the analysis}
\item{gene.names}{if "NULL", no gene name will be attached to the outputs,
otherwise it contains the vector of gene names}
\item{plot}{if "TRUE", plot the estimated pfp vs the rank of each gene}
\item{rand}{if specified, the random number generator will
be put in a reproducible state}
\item{calculateProduct}{if calculateProduct="TRUE"
(default) the rank product method is performed.
Otherwise the rank sum method is performed}
\item{MinNumOfValidPairs}{a parameter that indicates the minimum number
of NAs accepted per each gene. If it is set to NA (default) the half of
the number of replicates is used}
\item{RandomPairs}{number of random pairs generated in the function,
if set to NA (default),
the odd integer closer to the square of the number of replicates is used}
\item{huge}{if "TRUE" not all the outputs are evaluated
in order to save space}
\item{fast}{if "FALSE" the exact p-values for the Rank Sum are evaluated for
any size of the dataset.
Otherwise (default), if the size of the dataset is too big, only the p-values
that can be computed in "tail.time" minutes (starting from the tail) are
evaluated with the exact method. The others are estimated with the Gaussian
approximation. If calculateProduct="TRUE" this parameter is ignored}
\item{tail.time}{the time (default 0.05 min) dedicated to evaluate the exact
p-values for the Rank Sum.
If calculateProduct="TRUE" this parameter is ignored}
}
\value{
A summary of the results obtained by the Rank Product (or Rank Sum) method.
\item{pfp}{Estimated percentage of false positive predictions (pfp), both
considering upregulated an downregulated genes}
\item{pval}{Estimated pvalues per each gene being up- and down-regulated}
\item{RPs/RSs}{The rank-product/rank-sum statistics evaluated per each gene}
\item{RPrank/RSrank}{rank of the Rank Product (or Rank Sum) of each gene
in ascending order}
\item{Orirank}{Ranks obtained when considering each possible pairing.
In this version of the package, this is not used to compute
Rank Product (or Rank Sum), but it is kept for backward compatibility}
\item{AveFC}{Fold change of average expressions (class1/class2).
log fold-change if data has been log transformed,
original fold change otherwise}
\item{allrank1}{Fold change of class 1/class 2 under each origin.
log fold-change if data has been log transformed,
original fold change otherwise}
\item{allrank2}{Fold change of class 2/class 1 under each origin.
log fold-change if data has been log transformed,
original fold change otherwise}
\item{nrep}{Total number of replicates}
\item{groups}{Vector of labels (as cl)}
\item{RandomPairs_ranks}{a matrix containing the
ranks evaluated for each RandomPair}
}
\references{
Breitling, R., Armengaud, P., Amtmann, A., and Herzyk, P.(2004) Rank Products:
A simple, yet powerful, new method to detect differentially regulated genes in
replicated microarray experiments, FEBS Letter, 57383-92 
}
\author{
Francesco Del Carratore,
\email{francesco.delcarratore@postgrad.manchester.ac.uk}
\cr Andris Janckevics, \email{andris.jankevics@gmail.com}
}
\seealso{
\code{\link{topGene}} \code{\link{RP}} \code{\link{RPadvance}}
\code{\link{plotRP}} \code{\link{RP.advance}} \code{\link{RSadvance}}
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

RP.out <- RankProducts(golub[,subset],golub.cl[subset],rand=123) 
      
# class 2: label =1, class 1: label = 0
#pfp for identifying genes that are up-regulated in class 2 
#pfp for identifying genes that are down-regulated in class 2 
head(RP.out$pfp)
#Rank Sum
RS.out <-RankProducts(golub[,subset],golub.cl[subset],rand=123,
    calculateProduct=FALSE)
head(RS.out$pfp)
}
