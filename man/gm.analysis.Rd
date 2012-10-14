\name{gm.analysis}
\alias{gm.analysis}
\title{ Graphical model analysis }
\description{
  Graphical model selection by different selection strategies. 
  The selected graph is plotted where line thickness corresponds to an optional edge measure.
}
\usage{
gm.analysis(data, edge.measure = c("gamma.cond", "marg.gamma", "cor", "boot",
            "cov", "p.value"), program = c("gms","coco", "mim"),
            strategy = c("backwards", "forwards", "eh"), 
            plot.significant = TRUE, boot.N = 100, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Data frame or array. Variables need to be discrete and should have names. }
  \item{edge.measure}{ One of "gamma.cond", "marg.gamma", "cor", "boot", "cov", "p.value". To be plotted in the graph. May be abbreviated. See details. }
  \item{program}{ Selects whether the gamma selection \code{\link{gm.gms}}, MIM or the CoCo package is used for model selection. Calls either the function
                \code{\link{gm.coco}} or \code{\link{gm.mim}}. May be abbreviated. }
  \item{strategy}{ See \code{\link{gm.gms}}, \code{\link{gm.coco}} or \code{\link{gm.mim}}. }
  \item{plot.significant}{ If TRUE only significant edges in the selected models are plotted (in solid lines). 
                           If FALSE also not significant edges are plotted as dashed lines.
 }
  \item{boot.N}{ Only required if \code{edge.measure} is "boot". Specifies the number of bootstrap replications. }
  \item{\dots}{ Further options of \code{\link{gm.coco}} or \code{\link{gm.mim}} can be given. }
}
\details{
  Different edge measures:\cr
  "gamma.cond" - Conditional gamma coefficient, calls \code{\link{gm.gamma}} with type = conditional.\cr
  "marg.gamma" - Marginal gamma coefficient, calls \code{\link{gm.gamma}} with type = marginal.\cr
  "cor" - Pearson's correlation coefficient, calls \code{\link{cor}}.\cr
  "boot" - calls \code{\link{gm.boot.mim}} or \code{\link{gm.boot.coco}} depending on the choice of \code{program}. 
           \code{boot.N} bootstrap replications will be done. Plots the relative edge frequency.\cr
  "cov" - Pearson's covariance, calls \code{\link{cov}}. \cr
  "p.value" - is only available if \code{program} choice was "mim". Then \code{plot.significant} is set to TRUE since 
              the p-value is only available for edges present in the model. In this case the edges will be plotted thinner the smaller the p-value.
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
\seealso{ \code{\link{gm.plot}}, \code{\link{gm.coco}}, \code{\link{gm.mim}} }
\examples{
\dontrun{ 
## until CoCo is available again ...
  data(wam)
  gm.analysis(wam)
  
  ### showing various options in action
  gm.analysis(wam,program="c",strategy="f",edge.measure="b",
    boot.N=50,plot.significant=FALSE,recursive=TRUE,follow=TRUE,decomposable.mode=TRUE)
             
  gm.analysis(wam,program="m",edge.measure="p",options="u")
  }
}
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ graphs }
\keyword{ hplot }
\keyword{ models }
