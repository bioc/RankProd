\name{plotRP}
\alias{plotRP}
\title{Graphical Display of the Rank Product/Sum analysis}
\description{  
Plot a graph of the estimated pfp vs the 
number of identified genes}
\usage{
plotRP(x, cutoff=NULL)}
\arguments{
\item{x}{the value returned by function RP, 
RPadvance, RSadvance, RankProducts or RP.advance}
\item{cutoff}{The pfp threshold value used to select genes}
}  
\value{
A graphical display of the estimated pfp vs number of 
identified genes, which is also the gene rank of its 
original rank product/sum across all comparison. If cutoff 
is sepcified, a horizontal line will be plotted on the 
graphic to indicate the positon of the cutoff point, and 
all genes identified will be marked in red.
     
Two plots will be displayed, one for the identification of 
up-regulated genes in class 2, one for the identification of 
down-regulated genes in class 2
}
\author{Fangxin Hong, \email{fhong@salk.edu}
\cr Francesco Del
Carratore, \email{francesco.delcarratore@postgrad.manchester.ac.uk}}
\seealso{
\code{\link{topGene}}   \code{\link{RP}}  
\code{\link{RPadvance}} \code{\link{RSadvance}}
}
\examples{
# Load the data of Golub et al. (1999). data(golub) 
#contains a 3051x38 gene expression
# matrix called golub, a vector of length called golub.cl 
#that consists of the 38 class labels,
# and a matrix called golub.gnames whose third column contains the gene names.
data(golub)
 
#use a subset of data as example, apply the rank product method
subset <- c(1:4,28:30)
#Setting rand=123, to make the results reproducible,
#identify genes that are up-regulated in class 2 
#(class label =1)
RP.out <- RP(golub[,subset],golub.cl[subset], rand=123)
      
#plot the results
plotRP(RP.out,cutoff=0.05)
}