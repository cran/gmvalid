\name{gm.si}
\alias{gm.si}
\title{ Synergy Index}
\description{
  Computes the synergy index of two discrete variables in relation to 
a binary group variable.
}
\usage{
gm.si(X,Y,group,data=0,reference=c(1,1,2),conf.level = 0.95)
}
\arguments{
  \item{X}{ Index of the variable's position in \code{data} or a vector. }
  \item{Y}{ Index of the variable's position in \code{data} or a vector. }
  \item{group}{ Binary group or outcome variable addressed as variable index in 
                \code{data} or as vector.}
  \item{data}{ Data frame or table. }
  \item{reference}{ Vector to define the \code{reference} categories of \code{X}, \code{Y} and \code{group}.
                     By default, 2 is the reference category for \code{group}.}
  \item{conf.level}{ Confidence level of the interval (default is 0.95). }
}
\details{
    Rothman's synergy index (S) is an interaction measure between two
    discrete variables on a dichotomous outcome. 
The index equals 1 under additivity, S > 1 in the case of synergy and 
S < 1 in the presence of antagonism.
       
    The synergy index is originally constructed on (2 x 2)-tables, but (i x j)-tables can
    also be analyzed. Then (i-1) x (j-1) synergy indices are computed and combined to 
an overall synergy index.    
    
    The confidence intervalls are calculated using the asymptotic variance given in 
    Rothman (1974). 
}
\value{
A list containing:
  \item{OddsRatio }{A odds ratio table. }
  \item{covariance}{Covariance matrix for category specific synergy indices.
                    Not written if \code{X} or \code{Y} are binary.}
  \item{measure }{Matrix containing the estimate(s), standard deviation(s), 
                  confidence interval(s) and p-value(s). Figures in brackets show the
                  reference category respectively the category under consideration.\cr                 
                  If both factors \code{X} and \code{Y} are binary, confidence 
                  intervals for case-control as well as cohort designs are computed.
                  If at least one factor has more than 2 categories, the overall synergy index
                  with its corresponding confidence interval is computed that follows a 
                  case-control design.}
}
\references{ 
  Rothman K (1974) 
  \emph{The estimation of synergy or antagonism.}
  American Journal of Epidemiology, 103(5):506-511

  Rothman K (1986) 
  \emph{Modern Epidemiology.}
  Little, Brown and Company, Boston/Toronto.
 }
\author{ 
  Ronja Foraita, Fabian Sobotka \cr
  Bremen Institute for Prevention Research and Social Medicine \cr
  (BIPS)  \url{http://www.bips.uni-bremen.de}
   }
\note{ 
  It can occur that certain combinations of categories lead to a negative synergy indices.
  In that case confidence intervals cannot be computed. If so, use the \code{reference} option
  to re-order the categories of the variable(s) in question (see example below).
} 
\seealso{ \code{\link{gm.csi}} }
\examples{
  data(wynder)
  gm.si(1,2,3,wynder)

  # Smoking and alcohol in relation to oral cancer among male veterans under age 60.
  # (from "Modern Epidemiology")
  oral <- array(c(20,3,18,8,12,6,166,225),dim=c(2,2,2), 
            dimnames=list(Group=c("control","cases"),
            Smoker=c("no","yes"),Alcohol=c("no","yes")))
  oral.df <- expand.table(oral)
  # grouping variable is first in data frame
  gm.si(2,3,1,oral.df)
  
  # Effects must be ascending in respect to the reference category
  show.effect <- array(c(1,7,2,7,7,12,106,48),dim=c(2,2,2),
                        dimnames=list(A=1:2,B=1:2,C=1:2))
  # produces NaN
  gm.si(1,2,3,expand.table(show.effect))
  # > re-ordering variable B helps
  gm.si(1,2,3,expand.table(show.effect),reference=c(1,2,2))

}
\keyword{ htest }
