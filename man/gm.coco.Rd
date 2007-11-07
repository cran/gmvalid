\name{gm.coco}
\alias{gm.coco}
\title{ Graphical model selection (CoCo) }
\description{
  Selects one or more graphical models (depending on strategy) using the CoCo package.
  Different types of model selection are available. A base model can be optionally defined.
}
\usage{
gm.coco(data, strategy = c("backwards", "forwards", "eh"), model = FALSE, ...)
}
\arguments{
  \item{data}{ Data frame or a table (array). Variables should have names, \code{data} has to be discrete. }
  \item{strategy}{ Type of model selection. "backwards" searches for not significant edges to delete, starting from the saturated model as default (\code{\link[CoCoRaw]{backward}}).
                    "forwards" adds significant edges, starting from the main effects model (\code{\link[CoCoRaw]{forward}}). The "Edwards-Havranek"
                    model search rejects complete models in every step and finishes with one or more accepted models (\code{\link[CoCoRaw]{eh}}).
                    The default strategy is "backwards". Selections may be abbreviated.
                     }
  \item{model}{ Character string specifying a start model for "backwards" and "forwards" selection procedure.
                For "eh" a character string would specify the base model for tests.
                Since "eh" searches for models between a defined minimal and a maximal model, a string vector of length 2 is another option.
                In this case a minimal and a maximal model is to be given.
                The model formulae has to start with the first lowercase letters of the alphabet, e.g. "abc,cde". Variable names cannot be given.
                 }
  \item{\dots}{ Additional options to be given to the selection procedures. See \code{\link[CoCoRaw]{backward}} (forward) or \code{\link[CoCoRaw]{eh}}. }
}
\value{
A list containing:
  \item{"accepted" }{ String vector of the accepted model(s). }
  \item{"rejected" }{ String vector of the rejected model(s). }
  \item{"base"}{ Character string giving the base model the selection started with. }
  \item{"strategy"}{ Character string indicating what type of selection strategy was performed. }
  \item{"tests"}{ Matrix containing the test results of the accepted vs. the rejected models 
  (available if at least one of the models is submodel of another one).
                    Otherwise the accepted and rejected models are tested against the saturated model. }
  \item{"variable names"}{ Matrix that assigns a letter to each variable used in the model formulas. }
}
\author{   
  Ronja Foraita, Fabian Sobotka \cr
  Bremen Institute for Prevention Research and Social Medicine \cr
  (BIPS)  \url{http://www.bips.uni-bremen.de}
   }
\references{   Badsberg JH (2001)
            \emph{A guide to CoCo.}
            Journal of Statistical Software, 6(4).
             }
\note{ CoCo may crash R occasionally on a bad day. We tried to prevent all crash causes for discrete data. The R version of CoCo can only handle one CoCoObject.
        You may run "endCoCo()" if you use CoCo also outside this package or the function aborts.
}
\seealso{ \code{\link{gm.mim}}, \code{\link[CoCoRaw]{backward}}, \code{\link[CoCoRaw]{eh}} }
\examples{
  data(wam)
  gm.coco(wam)
  ### giving many options to the strategy
  gm.coco(wam,recursive=TRUE,follow=TRUE,decomposable.mode=TRUE,
            coherent=TRUE,IC=TRUE,BIC=TRUE)
  
  ### giving base
  gm.coco(wam,strategy="e",model=c("ab,cd","ae,be"))
  gm.coco(wam,strategy="f",model="abc,cd,de,f")
}
\keyword{ models }
\keyword{ multivariate }
