\name{arab}
\docType{data}
\alias{arab}
\alias{arab.cl}
\alias{arab.origin}
\alias{arab.gnames}
\title{Genomic Response to Brassinosteroid in Arabidopsis}
\description{
       This data if from Affy ATH1 array experiments of 
       genomic response to brassinosteroid in Arabidopsis 
       conducted at two laoratories. The data set contains 
       random selected 500 genes and 10 samples, 6 from lab 1 and 4 from 
       lab 2. Data was pre-processed by RMA 
}

\usage{
data{arab}
}

\format{ 
   \item{arab}{matrix of gene expression levels of 500 genes 
               from 10 samples, rows correspond to genes  and 
               columns to mRNA samples.}
   \item{arab.cl}{numeric vector indicating the treatment class, 
                 5 brassinosteroid-treated cases (code 1) 
                  and 5 control cases (code 0) }
   \item{arab.gnames}{character vector containing the AffyID of 
                  the 500 genes for the expression matrix \code{arab}}
   \item{arab.origin}{numeric vector indicating the origin of 
                     the samples, 6 samples from lab 1 (code 1) 
                     and 4 samples from lab 2 (code 2)} 
} 
\references{
    Nemhauser JL, Mockler TC, Chory J.  Interdependency of 
    brassinosteroid and auxin signaling in Arabidopsis. 
    \emph{PLoS Biol.} 2004 21460.
  
    Microarray data from AtGenExpress 
    (http://arabidopsis.org/info/expression/ATGenExpress.jsp)

}
\keyword{datasets}    
        