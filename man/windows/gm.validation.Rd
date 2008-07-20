\name{gm.validation}
\alias{gm.validation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Validation and variability measures for graphical models. }
\description{
  The bootstrapped graphical models are analyzed and some new uncertainty measures
  and variance estimators are applied in order to determine the uncertainty of a
  selected graphical model.
}
\usage{
gm.validation(data, N = 0, program = c("coco", "mim"), Umax = 0.5, conf.level = 0.95, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Either the output list from one of the bootstrap functions \code{gm.boot.coco} or \code{gm.boot.mim}
               or data frame or a table (array). Variables should have names, \code{data} has to be discrete. }
  \item{N}{ Number of bootstrap replications. Only needed if data is not yet a bootstrap output. }
  \item{program}{ If a bootstrap is not yet done, which function should do it: \code{gm.boot.coco} or \code{gm.boot.mim}. }
  \item{Umax}{ Defines the maximum uncertainty in the edges. In the validation the selection frequency of every edge in the
               bootstrap replications is the most important foundation. The default is that an edge is maximally uncertain
               when it is selected in 50% of the replications.
                }
  \item{conf.level}{ Confidence level for the bootstrap percentile interval.
                    }
  \item{\dots}{ You may add options to the selection strategy, if a bootstrap is still needed. 
                See \code{\link{gm.boot.coco}} or \code{\link{gm.boot.mim}}. 
                }
}
\details{
  The bootstrap functions bring multivariate output about the uncertainty of a selected graphical model.
  This function presents some possibilities to reduce the uncertainty to a univariate measure, based either on
  the edge frequencies of presence in the bootstrapped models or on differences between models measured in edges.
}
\value{
  \item{ "original model" }{ Character string giving the selected model using the original unsampled data. }
  \item{ "mode model" }{ Character string giving a model that was selected the most as a whole.
                              It can happen that for 100 bootstrap replications there are 100 different models selected.
                              So this is not the best way to validate a graphical model.
                               }
  \item{ "mean model" }{ Character string giving a model that consists of all the edges that had a frequency of selection 
                              greater than \code{Umax} in the bootstrap replications. This selection promisses the best
                              validation results.
                              }
  \item{ "MEU" }{ Linear measure for the uncertainty of the \code{validated model} based on the edge frequency.
                                  For each edge the minimum distance to a frequency of 0 or 1 is measured, multiplied by \code{Umax}
                                  and identified as its uncertainty because with a frequency of 0 or 1 we think of an edge as maximally safe.
                                  Then the mean over all edges is calculated to receive an uncertainty measure.
                                   }
  \item{ "MSEU" }{ Same as above, only the uncertainty of an edge is squared and normalized
                                          so that a frequency of \code{Umax} results in a stability of 0 respectively
                                          an uncertainty of 1.
                                           }
  \item{ "differing edges" }{ List that says by how many edges the bootstrapped models differ from the \code{validated model}. }
  \item{ "total possible edges" }{ Number of edges the saturated model would include. }
  \item{ "model std" }{ The validated model is interpreted as mean model and the number of differing edges
                             define the difference between 2 models and so the univariate standard error of the \code{validated model}
                             is estimated. }
  \item{ "std/total" }{ Uncertainty measure that relates the standard error (\code{model std}) to the value of \code{total possible edges}. }
  \item{ "expected edges different" }{ The calculated mean of the \code{differing edges}. }
  \item{ "bootstrap percentile 95" }{ The value that includes at least the lower 95% of \code{differing edges} to give an upper border to
                                        the model uncertainty.
                                        }
  \item{ "variable names" }{ Matrix that assigns a letter to each variable that is used in the model formulas. }
}
\references{ 
            Efron B, Tibshirani RJ (1993) 
            \emph{An Introduction to the Bootstrap.}
             Chapman & Hall
              }
\author{ 
  Ronja Foraita, Fabian Sobotka \cr
  Bremen Institute for Prevention Research and Social Medicine \cr
  (BIPS)  \url{http://www.bips.uni-bremen.de}
   }
\note{ 
  The question when an edge is maximally uncertain is not yet answered satisfactory.
  Can we say that an association that is selected in 40% of the bootstrap replications
  is not present? Or that an edge that is present in 60% of the cases is not significant
  randomly? Therefore the argument \code{Umax} leaves it up to your opinion.
  
  If you already have run a bootstrap, make sure it was with all the possible \code{calculations}.
}
\seealso{ \code{\link{gm.boot.coco}}, \code{\link{gm.boot.mim}} }
\examples{
  ### Standard procedure
  data(wam)
  boot.out <- gm.boot.coco(1000,wam,strategy="f",recursive=TRUE,follow=TRUE,all.significant=FALSE)
  gm.validation(boot.out)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ models }% __ONLY ONE__ keyword per line
\keyword{ nonparametric }
