\name{plotGene}
\alias{plotGene}
\title{Graphical Display of The expression levels}
\description{  
     Plot the expression values and print statistical results 
     for each individual gene based on a user query}

\usage{
     plotGene(gene.to.plot,x,gene.names=NULL, data,cl,origin,
              logged=TRUE,logbase=2)}

\arguments{
     \item{Gene.to.plot}{Name of the gene to be plotted}
     \item{x}{the value returned by function RP, 
           RPadvance or RSadvance }
     
     \item{gene.name}{Names of all genes in the data set. If "NULL", rownames 
             of the data will be used}
     \item{data.cl}{the same as that used in RP or RPadvance}

    \item{origin}{a vector containing the origin labels of the 
              sample. e.g. for 
            the data sets generated at multiple laboratories, the label
            is the same for samples within one lab and different for samples 
            from different labs. The same as that used in RPadvance} 
     \item{logged}{if "TRUE", data has bee logged, otherwise set it 
               to "FALSE"}

     \item{logbase}{base used when taking log, used to restore the 
                 fold change.The default value is 2, this will be 
                 ignored if logged=FALSE}


}  
     
\value{
     A graphical display of the expression levels of the input gene. 
     The estimated statistics for differential expression will be printed 
     on the plot as well as output in the screen.The statistics include: 
     F.C.:fold-change under each dataset if multiple datasets are used
     AveFC: average fold-change across all datasets
     pfp(pval): estimated percentage of false prediction (p-value) for 
       differential expression under each of the two tests: up-regulation 
       in class 2 compared with classs 1 and down-regulation in class 2 
       compared with class 1  
     
}

\author{Fangxin Hong \email{fhong@salk.edu}}

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
\keyword{htest}


