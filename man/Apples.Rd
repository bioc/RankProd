\name{Apples}
\docType{data}
\alias{apples.data}
\alias{apples.data.vsn}
\alias{apples.cl}
\alias{Biom}
\alias{mz}
\alias{rt}

\title{Metabolomics data on spiked apples}
\description{
A dataset of LC-MS features, obtained from twenty apples. The last ten apples
are spiked with known compounds. This set provides a test case for biomarker
selection methods: the task is to retrive the true biomarker variables.
The raw LC-MS data hava been converted to CDF format and processed with
XCMS to obtain the basepeaks.}

\usage{data(Apples)}

\value{
The format is a list of four elements:
\item{mz}{the m/z values of the features (rounded)}
\item{rt}{the retention times of the features}
\item{apples.data}{a matrix containing the intensities
in the individual samples}
\item{apples.data.vsn}{a matrix containing the intesities after variance
stabilization and normalization performed with the vsn package}
\item{Biom}{the indices of the "true" biomarkers}
\item{apples.cl}{numeric vector encoding which samples are part of the spiked
class (code 1) and which ones are controls (code 0)}
}
\author{
Francesco Del Carratore
}

\references{
P. Franceschi, D. Masuero, U. Vrhovsek, F. Mattivi and R. Wehrens: A benchmark
spike-in data set for biomarker identification in metabolomics.
\emph{J. Chemom}. 26, 16-24 (2012)

R. Wehrens, P. Franceschi, U. Vrhovsek and F. Mattivi. Stability-based biomarker
selection. \emph{Analytica Chimica Acta} (2011), 705, 15-23.
\url{http://dx.doi.org/10.1016/j.aca.2011.01.039}
}
\examples{
data(Apples)
## show features identified in all apples
plot(rt, mz,
     xlab = "Retention time (s)", ylab = "m/z",
     main = "Spiked apples - subset")
}

\keyword{datasets}    
        