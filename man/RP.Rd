\name{RP}
\alias{RP}
\title{Rank Product Analysis}
\description{
The function performs the Rank Product method to
identify differentially expressed genes. It is possible to do either a
one-class or two-class analysis. This function has been kept only to guarantee
backward compatibility, in fact the same results can be obtained by
\code{\link{RankProducts}}.}
\usage{
RP(data, cl, num.perm = 100, logged = TRUE, na.rm = TRUE, gene.names = NULL,
plot = FALSE, rand = NULL, huge = FALSE)
}
\arguments{
\item{data}{the function performs the Rank Product (or Rank Sum) method to
identify differentially expressed genes. It is possible to do either a
one-class or two-class analysis}
\item{cl}{a vector containing the class labels of the samples.
In the two class unpaired case, the label of a sample
is either 0 (e.g., control group) or 1 (e.g., case group).
For one class data, the label for each sample should be 1}
\item{num.perm}{in this version of the package, this parameter is not used any
more, but it is kept for backward compatibility}
\item{logged}{if "TRUE" data have been previously log transformed. Otherwise it 
should be set as "FALSE"}
\item{na.rm}{if "FALSE", the NA value will not be used in
computing rank. If "TRUE" (default), the missing values will be replaced by
the genewise median of the non-missing values.
Gene with a number of missing values greater than 50\% are still
not considered in the analysis}
\item{gene.names}{if "NULL", no gene name will be attached to the outputs,
otherwise it contains the vector of gene names}
\item{plot}{if "TRUE", plot the estimated pfp vs the rank of each gene}
\item{rand}{if specified, the random number generator will
be put in a reproducible state}
\item{huge}{if "TRUE" not all the outputs are evaluated
in order to save space}
}
\value{
A summary of the results obtained by the Rank Product method.
\item{pfp}{estimated percentage of false positive predictions (pfp), both
considering upregulated an downregulated genes}
\item{pval}{estimated pvalues per each gene being up- and down-regulated}
\item{RPs}{the Rank Product statistics evaluated per each gene}
\item{RPrank}{rank of the Rank Product of each gene
in ascending order}
\item{Orirank}{ranks obtained when considering each possible pairing.
In this version of the package, this is not used to compute
Rank Product (or Rank Sum), but it is kept for backward compatibility}
\item{AveFC}{fold changes of average expressions (class1/class2).
log fold-change if data has been log transformed,
original fold change otherwise}
\item{allrank1}{fold change of class 1/class 2 under each origin.
log fold-change if data has been log transformed,
original fold change otherwise}
\item{allrank2}{fold change of class 2/class 1 under each origin.
log fold-change if data has been log transformed,
original fold change otherwise}
\item{nrep}{total number of replicates}
\item{groups}{vector of labels (as cl)}
\item{RandomPairs_ranks}{a matrix containing the
ranks evaluated for each RandomPair}
}
\references{
Breitling, R., Armengaud, P., Amtmann, A., and Herzyk, P.(2004) Rank Products:A 
simple, yet powerful, new method to detect differentially regulated genes in 
replicated microarray experiments, FEBS Letter, 57383-92 
}
\author{
Francesco Del Carratore,
\email{francesco.delcarratore@postgrad.manchester.ac.uk}
\cr Andris Janckevics, \email{andris.jankevics@gmail.com}
}
\note{
Percentage of false prediction (pfp), in theory,
is equivalent of false discovery 
rate (FDR), and it is possible to be large than 1.

The function looks for up- and down- regulated genes in two seperate steps, thus 
two pfps and pvalues are computed and used to identify gene that belong to each 
group.

This function is suitable to deal with data from a single origin, e.g. single 
experiment. If the data has different origin, e.g. generated at different
laboratories, please refer RP.advance. 
}



\seealso{
\code{\link{topGene}} \code{\link{RSadvance}} \code{\link{RPadvance}}
\code{\link{plotRP}} \code{\link{RP.advance}} \code{\link{RankProducts}}
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