\name{gm.analysis}
\alias{gm.analysis}
\title{ Graphical model analysis }
\description{
  Graphical model selection by different selection strategies. 
  The selected graph is plotted where line thickness corresponds to an optional edge measure.
}
\usage{
gm.analysis(data, edge.measure = c("gamma.cond", "marg.gamma", "cor", "boot",
            "cov"), strategy = c("backwards", "forwards", "eh"), 
            plot.significant = TRUE, boot.N = 100, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Data frame or array. Variables need to be discrete and should have names. }
  \item{edge.measure}{ One of "gamma.cond", "marg.gamma", "cor", "boot", "cov". To be plotted in the graph. May be abbreviated. See details. }
  \item{strategy}{ See \code{\link{gm.coco}}. }
  \item{plot.significant}{ If TRUE only significant edges in the selected models are plotted (in solid lines). 
                           If FALSE also not significant edges are plotted as dashed lines.
 }
  \item{boot.N}{ Only required if \code{edge.measure} is "boot". Specifies the number of bootstrap replications. }
  \item{\dots}{ Further options of \code{\link{gm.coco}} can be given. }
}
\details{
  Different edge measures:\cr
  "gamma.cond" - Conditional gamma coefficient, calls \code{\link{gm.gamma}} with type = conditional.\cr
  "marg.gamma" - Marginal gamma coefficient, calls \code{\link{gm.gamma}} with type = marginal.\cr
  "cor" - Pearson's correlation coefficient, calls \code{\link{cor}}.\cr
  "boot" - calls \code{\link{gm.boot.coco}}. \code{boot.N} bootstrap replications will be done. 
           Plots the relative edge frequency.\cr
  "cov" - Pearson's covariance, calls \code{\link{cov}}. \cr
}
\value{
A list containing:
  \item{"strategy"}{ One of "backwards","forwards" and "eh". }
  \item{"model" }{ String vector of the accepted and plotted models. }
  \item{"edge.measure" }{ Character string indicating what type of edge measure was performed. }
  \item{"analysis"}{ Upper-tri matrix of the respective edge measures. }
  \item{"variable names"}{ Matrix that assigns a letter to each variable that is used in the model formulas. }
}
\author{ 
  Fabian Sobotka, Marc Suling, Ronja Foraita \cr
  Bremen Institute for Prevention Research and Social Medicine \cr
  (BIPS)  \url{http://www.bips.uni-bremen.de}
   }
\seealso{ \code{\link{gm.plot}}, \code{\link{gm.coco}} }
\examples{
  data(wam)
  gm.analysis(wam)
  
  ### showing various options in action
  gm.analysis(wam,strategy="f",edge.measure="b",
    boot.N=50,plot.significant=FALSE,recursive=TRUE,follow=TRUE,decomposable.mode=TRUE)
}
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ graphs }
\keyword{ hplot }
\keyword{ models }
