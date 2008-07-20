\name{gm.gms}
\alias{gm.gms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Graphical model selection for nominal / ordinal data using Goodman's and Kruskal's Gamma. }
\description{
  Selects a graphical model forwards or backwards in one or more steps using the
  conditional gamma coefficient impelemented in \code{\link{gm.gamma}}.
}
\usage{
gm.gms(data, strategy = c("backwards", "forwards", "combined"), model = FALSE, onestep = FALSE,
       headlong=FALSE, conf.level = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Data frame or table (array). Variables should have names. \code{data} has to be discrete. }
  \item{strategy}{ Type of model selection. "backwards" searches for not significant edges to delete, starting from the saturated model as default.
                    "forwards" adds significant edges, starting from the main effects model. 
                    The "combined" strategy is a 3 step procedure: \code{\link{gm.screening}}, "backwards" and then "forwards".
                    The default strategy is "backwards". Selections may be abbreviated.
                     }
  \item{model}{ Character string specifying a start model for "backwards" and "forwards" selection procedure.
                The model formulae has to start with the first lowercase letters of the alphabet, e.g. "abc,cde". Variable names cannot be given.
                 }
  \item{onestep}{ If TRUE all edges associated with a p-value < 1 - conf.level / p-value > 1 - conf.level will be added / removed
                  in one step. If FALSE only one edge will be added / removed in each step. }
  \item{headlong}{ If TRUE the edges are visited in random order and the first (in)significant is removed / added. Helps to avoid making the
                    easiest decisions about edges in the first steps and later not being able to make hard decisions and the model ends up having
                    too many edges. If FALSE in every step the edge with the highest / lowest p-value is removed / added.
                    Only working for \code{onestep} FALSE.
                    }
  \item{conf.level}{ See \code{\link{gm.gamma}}. }
}
\details{
  For every two-variable association the conditional gamma coefficient (in case of no other dependencies the marginal gamma),
  the standard deviation and the p-value is calculated from the data. Backwards the association with the smallest p-value is kicked
  out in every step. Forwards the highest is included until there are no more significances. In the one-step procedure all (in-)
  significant are added / deleted at once, where the basis is the marginal / saturated model while when doing more steps the base model
  is always the one from the previous step.
}
\value{
  \item{measure }{ A list of matrices with the output of gm.gamma for all cliques in the selected model. }
  \item{model }{ The selected model in a string. }
}
\references{
  Davis JA (1967) 
  \emph{A partial coefficient for Goodman and Kruskal's gamma.}
  Journal of the American Statistical Association, 62:189-193.

  Olszak M, Ritschard G (1995) 
  \emph{The behaviour of nominal and ordinal partial association measures.}
  The Statistician, 44(2):195-212.
   }
\author{
  Ronja Foraita, Fabian Sobotka \cr
  Bremen Institute for Prevention Research and Social Medicine \cr
  (BIPS)  \url{http://www.bips.uni-bremen.de}
   }
\note{ 
       The function is more time consuming than comparable functions.
}
\seealso{ \code{\link{gm.coco}}, \code{\link{gm.gamma}} }
\examples{

data(wam)
gm.gms(wam)

gm.gms(wam,onestep=TRUE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
\keyword{ nonparametric }
