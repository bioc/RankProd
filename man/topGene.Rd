
\name{topGene}
\alias{topGene}
\title{Output Significant Genes}
\description{  
     Identify differentially expressed genes using 
         rank product method}

\usage{

     topGene(x,cutoff=NULL,method="pfp",num.gene=NULL,logged=TRUE,logbase=2,gene.names=NULL)}

\arguments{

     \item{x}{the value returned by the 
           function RP, RP.advance or
          Rsum.advance }
   
     \item{cutoff}{threshold in pfp used to select genes}
  
      \item{method} {If cutoff is provided, the method needs to be 
              selected to identify genes."pfp" uses percentage of 
              false prediction, which is the default setting.
              "pval" used p-value which is less 
              stringent than pfp}

     \item{num.gene}{number of candidate genes of interests, 
            if cutoff is provided, this will be ignored}

     \item{logged}{if "TRUE", data has bee logged, otherwise set it 
               to "FALSE"}

     \item{logbase}{base used when taking log, used to restore the 
                 fold change.The default value is 2, this will be 
                 ignored if logged=FALSE}

     \item{gene.names}{if "NULL", no gene name will be 
                       attached to the output table} 

}


\value{

     Two tables of identified genes with 
     gene.index: index of gene in the original data set 
     RP/Rsum: Computed rank product/sum for each gene
     FC:(class1/class2): Expression Fold change of class 1/ class 2.                   
     pfp: estimated pfp for each gene if the gene is used as cutoff point
     P.value: estimated p-value for each gene 

     Table 1 list genes that are up-regulated under class 2, Table 1 ist 
     genes that are down-regulated under class 2, 
}


\author{Fangxin Hong \email{fhong@salk.edu}}

\seealso{

    \code{\link{plotRP}} \code{\link{RP}}  
     \code{\link{RPadvance}} \code{\link{RSadvance}}

}


\reference{

      Breitling, R., Armengaud, P., Amtmann, A., and Herzyk, 
      P.(2004) Rank Products: A simple, yet powerful, new method 
      to detect differentially regulated genes in
      replicated microarray experiments, \emph{FEBS Letter}, 57383-92
}
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

      #identify genes 
      RP.out <- RP(golub[,subset],golub.cl[subset],rand=123)  

      #get two lists of differentially expressed genes 
      #by setting FDR (false discivery rate) =0.05

      table=topGene(RP.out,cutoff=0.05,method="pfp",logged=TRUE,logbase=2,
                   gene.names=golub.gnames[,3])
      table$Table1
      table$Table2

      #using pvalue<0.05
      topGene(RP.out,cutoff=0.05,method="pval",logged=TRUE,logbase=2,
                   gene.names=golub.gnames[,3])

      #by selecting top 10 genes

      topGene(RP.out,num.gene=10,gene.names=golub.gnames[,3])

}



\keyword{htest}

