\name{gm.screening}
\alias{gm.screening}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Model Screening }
\description{
  Preliminary model search procedure for categorial data in three steps, using 
  Goodman's and Kruskal's gamma and chi-squared tests.
}
\usage{
gm.screening(data, conf.level = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data frame or a table (array). Variables should have names, \code{data} has to be discrete. }
  \item{conf.level}{ Confidence level for tests on independence. }
}
\details{
  The initial screening divides into three parts. In the first part edges are added to the main effects model, 
  if the marginal gamma coefficient or the chi-squared statistic are significant. In the second part an edge is 
  added between two variables, if at least one of the conditional gamma
  coefficient and the conditional chi-squared statistic is significant 
  for any possible variable in condition. In the third part all triples of variables that build
  a complete subgraph are tested for removal using the conditional gamma coefficient and conditional chi-squared statistic. 
  An edge is removed if the corresponding p-value is larger than 1 - \code{conf.level}. 
}
\value{
  \item{ mat }{ Adjacency matrix of the screened model. }
  \item{ model }{ String of the screened model. }
}
\references{ 
    Kreiner S and Edwards D (1983)
    \emph{ The analysis of contingency tables by graphical models }
    Biometrika, 70(3):553-565
    
    Kreiner S (2008)
    \emph{ DIGRAM }
    \url{http://staff.pubhealth.ku.dk/~skm/skm/} last accessed 30.06.2008

    Siersma V (2007)
    \emph{ Studies in the interactions between disease development and interventions }
    PhD thesis, Faculty of Health Sciences, University of Copenhagen.

 }
\author{   
  Ronja Foraita, Fabian Sobotka \cr
  Bremen Institute for Prevention Research and Social Medicine \cr
  (BIPS)  \url{http://www.bips.uni-bremen.de}
   }
\note{ 
  It is advised to select a model after screening based on a backwards and forwards selection
  as implemented in \code{\link{gm.coco}}.
}
\seealso{ \code{\link[CoCoRaw]{backward}}, \code{\link{gm.coco}}, \code{\link{gm.chi}}, \code{\link{gm.gamma}} }
\examples{
  data(wam)
  model <- gm.screening(wam)
  gm.coco(wam,recursive=TRUE,criterion="aic",model=model$model)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
