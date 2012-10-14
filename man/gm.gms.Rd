\name{gm.gms}
\alias{gm.gms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Graphical model selection for nominal / ordinal data using Goodman's and Kruskal's Gamma. }
\description{
  Selection strategy for a graphical model. Searches forwards or backwards in one or more steps using the
  conditional gamma coefficient impelemented in \code{\link{gm.gamma}}.
}
\usage{
gm.gms(data, strategy = c("backwards", "forwards", "combined"), model = FALSE, onestep = FALSE, headlong = FALSE, conf.level = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Data frame or array. Variables need to be discrete and should have names. }
  \item{strategy}{ Type of model selection. "backwards" eliminates not significant edges, starting from the saturated model as default.
                    "forwards" adds significant edges, starting from the main effects model.
                    The "combined" strategy is a 3 step procedure: \code{\link{gm.screening}}, "backwards" and then "forwards".
                    The default strategy is "backwards". Selections may be abbreviated.
                     }
  \item{model}{ Character string to specify a start model for "backwards" or "forwards" selection procedures.
                For "combined" the model cannot be given, a start model 
                will be specified by \code{\link{gm.screening}}.
                The model formula has to start with the first lowercase letters of the alphabet, e.g. "abc,cde". 
                Variable names cannot be given.
                 }
  \item{onestep}{ If TRUE all edges associated with a p-value < 1 - conf.level / p-value > 1 - conf.level will be added / removed
                  in one step. If FALSE only one edge will be added / removed in each step. }
  \item{headlong}{ If TRUE edges are visited in random order and the first (in)significant is added / eliminated.
                   If FALSE in every step the edge with the highest / lowest p-value is eliminated / added.
                    Only working for \code{onestep} FALSE.
                    }
  \item{conf.level}{ See \code{\link{gm.gamma}}. }
}
\details{
  For every two-variable association the conditional gamma coefficient,
  the standard error and the p-value is calculated from the data. In the one-step procedure all (in-)
  significant edges are added / deleted at once, where the basis is the main effects / saturated model while when doing more steps the base model
  is always the selected model from the previous step.
}
\value{
  \item{measure }{ A list of matrices with the output of gm.gamma for all cliques in the selected model. }
  \item{model }{ Srting of selected model. }
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
\seealso{ \code{\link{gm.coco}}, \code{\link{gm.mim}}, \code{\link{gm.gamma}} }
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
