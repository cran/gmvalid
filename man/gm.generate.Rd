\name{gm.generate}
\alias{gm.generate}
\title{ Random data frames of binary variables given fixed marginals }
\description{
  Generates a random data frame of independent binary variables given fixed marginals.
}
\usage{
gm.generate(N, p = c(0.5, 0.5, 0.5))
}
\arguments{
  \item{N}{ Number of observations, sample size. }
  \item{p}{ Probability vector. Each element assigns the probability to draw a "2". The user-defined number of
            elements determines the number of variables in the resulting data frame. }
}
\value{
    A data frame with number of rows equal to \code{N} and number of columns equal to the length of \code{p}.
}
\author{ 
  Ronja Foraita, Fabian Sobotka \cr
  Bremen Institute for Prevention Research and Social Medicine \cr
  (BIPS)  \url{http://www.bips.uni-bremen.de}
 }
\seealso{\code{\link{gm.modelsim}}, \code{\link{gm.sim.ixj}}, \code{\link{r2dtable}} }
\examples{
gm.generate(10,c(.5,.2,.2))
gm.generate(15,c(.5,.5,.5,.5,.5,.5))
}
\keyword{ datagen }
\keyword{ distribution }
