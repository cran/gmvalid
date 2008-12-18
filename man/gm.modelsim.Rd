\name{gm.modelsim}
\alias{gm.modelsim}
\title{ Random data frames with given dependence model and fixed marginals }
\description{
    Generates a random data frame of discrete variables given a dependence model and marginals.
}
\usage{
    gm.modelsim(N, model, categories = 0)
}
\arguments{
  \item{N}{ Number of observations, sample size. }
  \item{model}{ A character string assigning a dependence model expressed as clique structure. Each variable
                has to be expressed as a letter, e.g. "ABC,CDE".
                             }
  \item{categories}{ a list of weight vectors that assigns the weight of each catogory.
                    Number of list elements must equal the number of variables in \code{model}.
                    Default is "list(c(.5,.5),c(.5,.5),...)".
                    }
}
\value{
    A data frame with number of rows approximately equal to \code{N} and number of columns equal
    to the number of variables used in \code{model}. 
}
\author{ 
  Ronja Foraita, Fabian Sobotka \cr
  Bremen Institute for Prevention Research and Social Medicine \cr
  (BIPS)  \url{http://www.bips.uni-bremen.de}
 }
\note{ 
    Observed marginal probabilities reflect the given marginal probabilites only approximatively.
    Works best with large sample sizes, like e.g. N=10,000.
}
\seealso{ \code{\link{gm.generate}}, \code{\link{gm.sim.ixj}}, \code{\link{r2dtable}} }
\examples{
    gm.modelsim(100,"AB,AC")
    table( gm.modelsim(100,"a,b,c") )
    
    tmp.df <- gm.modelsim(10000,"abf,cd,cf,bdeg,bfg")
    
    # with given number of categories
    tmp.df <- gm.modelsim(1000,"AB,C",list(c(1,1,1),c(1,1),c(1,1,1)))

    # with given number of categories and marginals
    tmp.df <- gm.modelsim(1000,"ABC",list(c(0.3,0.3,0.4),c(0.6,.4),c(0.25,0.25,0.5)))
    table(tmp.df)

    \dontrun{tmp.df <- gm.modelsim(100,"ABC",list(3,2,3))# (number of categories will be 2 x 2 x 2 )
            gm.modelsim(100,"123")
            }
}
\keyword{ datagen }
\keyword{ graphs }
