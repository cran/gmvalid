\name{gm.screening}
\alias{gm.screening}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Model Screening }
\description{
  Preliminary model search procedure for categorical data in three steps, using 
  Goodman's and Kruskal's gamma and chi-squared tests.
}
\usage{
gm.screening(data, conf.level = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data frame or array. Variables need to be discrete and should have names. }
  \item{conf.level}{ Confidence level (default is 0.95). }
}
\details{
  The initial model screening divides into three parts:
  (1) An edge is added to the main effects model, 
  if the marginal gamma coefficient or the chi-squared statistic is significant. 
  (2) An edge is added between two vertices, if either the conditional gamma
  coefficient or the conditional chi-squared statistic is significant 
  for the two variables given any possible variable in condition. 
  (3) All triplets of variables whose vertices build
  a complete subgraph are tested for removal using the conditional gamma coefficient and conditional chi-squared statistic. 
  An edge is removed if the corresponding p-value is larger than 1 - \code{conf.level}. 
}
\value{
  \item{ mat }{ Adjacency matrix of the screened model. }
  \item{ model }{ String of the graphical model. }
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
  The model screening is not a model selection strategy! It 'just' searches for a good 
  start model to initialize a model selection. Common usage is to select 
  a model after screening based on a backwards and forwards selection
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
