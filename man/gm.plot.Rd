\name{gm.plot}
\alias{gm.plot}
\title{ Plot graphical models }
\description{
  Plots given graphical models and writes provided measures next to the edges.
}
\usage{
gm.plot(model, significant = TRUE, data.analysis = 0)
}
\arguments{
  \item{model}{ String vector with model formulas. See \code{\link{gm.modelsim}}. }
  \item{significant}{ If TRUE only significant edges in the selected models are plotted (in solid lines). If FALSE also not significant edges are plotted as dashed lines.
   }
  \item{data.analysis}{ Upper-tri matrix with the measure for the edge between variables i and j (i>j) at matrix position [i,j]. If the length of the \code{model} is bigger than 1, \code{data.analysis} has to be a list of matrices. }
}
\details{
  The line width of the edges will depend on the size of the numbers in \code{data.analysis}.
}
\value{
  TRUE
}
\author{ 
  Fabian Sobotka, Marc Suling, Ronja Foraita \cr
  Bremen Institute for Prevention Research and Social Medicine \cr
  (BIPS)  \url{http://www.bips.uni-bremen.de}
 }
\note{ Every use of the plot function opens a new window.
}
\seealso{ \code{\link{gm.analysis}} }
\examples{
  gm.plot("ABC,CDE")

  gm.plot("VBA,EVC")
  
  gm.plot(c("ABC,CDE","AB,BC,CD,DE","ABC,DEF,GHI"))
  
  gm.plot("AB,AC",FALSE,matrix(0.5,nrow=3,ncol=3))
}
\keyword{ hplot }
\keyword{ graphs }
