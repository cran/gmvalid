\name{gm.sim.ixj}
\alias{gm.sim.ixj}
\title{ Random (i x j)-way dependency tables with given marginals }
\description{
    Generates random (i x j)-way dependency tables with given marginals.
}
\usage{
gm.sim.ixj(N, pa, pb)
}
\arguments{
  \item{N}{ Number of observations; sample size. }
  \item{pa}{ Weight vector of length i to assign the distribution of the row totals. }
  \item{pb}{ Weight vector of length j to assign the distribution of the column totals. }
}
\value{
    A (i x j)-matrix.
}
\author{ 
  Ronja Foraita, Fabian Sobotka \cr
  Bremen Institute for Prevention Research and Social Medicine \cr
  (BIPS)  \url{http://www.bips.uni-bremen.de}
   }
\note{     
    Observed totals that were generated according to the shorter weight vector 
    reflect the given weights only approximatively.
}
\seealso{ \code{\link{gm.generate}}, \code{\link{gm.modelsim}}, \code{\link{r2dtable}} }
\examples{
    gm.sim.ixj(1000,c(1,1,1),c(1,1,1,1,1))
    gm.sim.ixj(1000,c(.2,.3,.4,.1),c(.5,.2,.3)) 
    
    round(gm.sim.ixj(30,c(1,1),c(1,1)))       
    
    tab <- round(gm.sim.ixj(500,c(.5,.5),c(.5,.5)))
    chisq.test(tab)   
}
\keyword{ datagen }
