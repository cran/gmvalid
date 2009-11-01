\name{gm.validation}
\alias{gm.validation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Validation and uncertainty measures for graphical models. }
\description{
  The bootstrapped graphical models are analyzed and some new uncertainty measures
  are applied in order to determine the uncertainty of a
  selected graphical model.
}
\usage{
gm.validation(data, N = 0, Umax = 0.5, conf.level = 0.95, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Output list from \code{gm.boot.coco} or \code{gm.boot.mim},
               data frame or array. Variables need to be discrete and should have names. }
  \item{N}{ Number of bootstrap replications. Only needed if data is not yet a bootstrap output. }
  \item{Umax}{ Parameter that defines the maximum uncertainty in the edge selection frequency (default is 0.5).
                }
  \item{conf.level}{ Confidence level for bootstrap percentile interval (default is 0.95).
                    }
  \item{\dots}{ To add options for the model selection strategy. Only needed if data is not yet a bootstrap output.
                See \code{\link{gm.boot.coco}}. 
                }
}
\details{
  The bootstrap functions bring multivariate output about the uncertainty of a selected graphical model.
  This function presents some possibilities to reduce the uncertainty to a univariate measure, based either on
  the edge frequencies of presence in the bootstrapped models or on differences between models measured in edges.
}
\value{
  \item{ "original model" }{ Character string of the selected graphical model, using the original unsampled data. }
  \item{ "mode model" }{ Character string of the model that was selected most frequently.
                               }
  \item{ "mean model" }{ Character string of a graphical model that consists of those edges whose selection frequency is 
                              greater than \code{Umax} over the bootstrap replications.
                              }
  \item{ "MEU" }{ Mean edge uncertainty. Linear measure for the uncertainty of the \code{mean model} based on the edge frequency \eqn{f^*(e)}{f*(e)}.
                  \deqn{MEU(\bar{G}^*) = \frac{1}{|V \times V|} \sum_{e \in V \times V}^{}
                  {\frac{\min\left[\frac{1-U_{max}}{0.5}f^*(e),\frac{U_{max}}{0.5}\left(1-f^*(e)\right)\right]}
                  {\frac{U_{max}(1-U_{max})}{0.5^2}}}}{%
                  MEU = 1/|VxV| sum[min[(1-Umax)/0.5 f*(e),Umax/0.5 (1-f*(e))] / [Umax (1-Umax) / 0.5 / 0.5 ]]}
                  with vertices set V.
                                   }
  \item{ "MSEU" }{ Mean squared edge uncertainty.
                   \deqn{MSEU(\bar{G}^*) = \frac{1}{|V \times V|} \sum_{e \in V \times V}^{}
                   {\left(\frac{\min\left[\frac{1-U_{max}}{0.5}f^*(e),\frac{U_{max}}{0.5}\left(1-f^*(e)\right)\right]}
                   {\frac{U_{max}(1-U_{max})}{0.5}}\right)^2}}{%
                   MSEU = 1/|VxV| sum[min[(1-Umax)/0.5 f*(e),Umax/0.5 (1-f*(e))] / [Umax (1-Umax) / 0.5  ]]^2}
                                           }
  \item{ "edge differences" }{ Frequency list of edges that differ in the bootstrapped models from the \code{mean model}. }
  \item{ "total possible edges" }{ Number of edges in the saturated model. }
  \item{ "model std" }{ \deqn{\widehat{std}_{boot}(\bar{G}^*) = \sqrt{\frac{1}{B-1} \sum_{b=1}^{B}{d(G_b^* , \bar{G}^*)^2} }}{%
                        std = 1/(B-1) sum[d(G*b,G*)^2]}
                        with number of bootstrap replications B, edge differences \eqn{d(G^*,\bar{G}^*)}{d(G*b,G*)}.  
                        }
  \item{ "MED" }{ Mean edge deviation. Mean of \code{edge differences}. }
  \item{ "bootstrap percentile 95" }{ The value that includes at least the lower 95\% of \code{edge differences} to give an upper border to
                                        the model uncertainty.
                                        }
  \item{ "variable names" }{ Matrix. Assigns a letter to each variable used in the model formulas. }
}
\references{ 
            Efron B, Tibshirani RJ (1993) 
            \emph{An Introduction to the Bootstrap.}
             Chapman & Hall
              
              
            Foraita R, Sobotka F, Pigeot I (2009)
            \emph{The uncertainty of a selected graphical model.}
            unpublished
            
            Sobotka F, Foraita R, Eberle A, Pigeot I (2008)
            \emph{GMVALID: An R-package to validate graphical models using data of external causes
                  of morbidity and mortality in Bremen, 1999-2006 (in german, poster)}
                  \url{http://www.bips.uni-bremen.de/ma_downloads/bips_gmvalid.pdf}
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
\seealso{ \code{\link{gm.boot.coco}} }
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
