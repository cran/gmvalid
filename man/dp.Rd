\name{dp}
\alias{dp}
\docType{data}
\title{ Death penalty example of Simpson's paradox }
\description{
  Well known example to illustrate Simpson's paradox.
  The data set shows that ignoring victim's race lead to a different conclusion
  than including victim's race in the analysis.
}
\usage{data(dp)}
\format{
  A data frame with 326 observations on the following 3 variables.
  \describe{
    \item{\code{Defendants.Race}}{A factor with levels \code{white} and \code{black}}
    \item{\code{Victims.Race}}{A factor with levels \code{white} and \code{black}}
    \item{\code{Death.Penalty}}{A factor with levels \code{yes} and \code{no}}
  }
}
\source{
Radelet ML (1981) \emph{Racial characteristics and the imposition of the Death penalty.}
American Sociological Review, 46(6):918-927.
}
\examples{
data(dp)
## Graphical model analysis shows that 'defendant's race' is 
## independent from 'death penalty' given 'victim's race'.
\dontrun{
## until CoCo is available again ...
gm.analysis(dp,program="coco",recursive=TRUE)
}
}
\keyword{datasets}
