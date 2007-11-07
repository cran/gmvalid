\name{idd35}
\alias{idd35}
\docType{data}
\title{ Type 1 Diabetes susceptibility loci Idd3 and Idd5 }
\description{
  Loci Idd3 and Idd5 are suspectable to be associated with type 1 diabetes.
  This data reflects an animal experiment from Hill et al. (2000).
}
\usage{data(idd35)}
\format{
  A data frame with 313 observations on the following 3 variables.
  \describe{
    \item{\code{idd3}}{a factor with levels \code{NN} and \code{BB}}
    \item{\code{idd5}}{a factor with levels \code{NN} and \code{BB}}
    \item{\code{group}}{a factor with levels \code{controls} \code{cases}}
  }
}
\details{
 Data was used from Cordell et al. (2001) to model the joint effect of alleles at
 the loci Idd3 and Idd5 on the outcome "type 1 diabetes". Only homozygeous genotypes 
 were available.
}
\source{
  Hill NJ, Lyons PA, Armitage N, Todd JA, Wicker LS and Peterson LB (2000)
  \emph{NOD Idd5 locus controls insulitis and diabetes and overlaps the 
  orthologous CTLA4/IDDM12 and NRAMP1 loci in humans.}
  Diabetes, 49:1744-1747.
}
\references{
  Cordell HJ, Todd JA, Hill NJ, Lord CJ, Lyons PA, Peterson LB, Wicker LS and Clayton DG (2001)
  \emph{Statistical modeling of interlocus interactions in a complex disease: Rejection
  of the multiplicative model of epistasis in type 1 diabetes.}
  Genetics, 158:357-367.
}
\examples{
  data(idd35)
  table(idd35)
}
\keyword{datasets}
