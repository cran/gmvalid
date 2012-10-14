\name{gm.coco}
\alias{gm.coco}
\title{ Graphical model selection (CoCo) }
\description{
  Selects one or more graphical models (depending on strategy) using the CoCo package.
  Different types of model selection are available. A base model can be optionally defined.
}
\usage{
gm.coco(data, strategy = c("backwards", "forwards", "eh", "combined"), 
        model = FALSE, eh.strategy = c("smallest", "alternating", "rough"), 
        criterion = c("lr","aic","bic"),tests=FALSE, ...)
}
\arguments{
  \item{data}{ Data frame or array. Variables need to be discrete and should have names. }
  \item{strategy}{ Type of model selection. "backwards" eliminates not significant edges, starting from the saturated model as default 
                    (\code{backward} from \pkg{CoCo}).
                    "forwards" adds significant edges, starting from the main effects model (\code{forward} from \pkg{CoCo}). The "eh" ("Edwards-Havranek")
                    model search rejects complete models in every step and finishes with one or more accepted models (\code{eh} from \pkg{CoCo}).
                    The "combined" strategy is a 3 step procedure: \code{\link{gm.screening}}, "backwards" and then "forwards".
                    The default strategy is "backwards". Selections may be abbreviated.
                     }
  \item{model}{ Character string specifying a start model for "backwards" or "forwards".
                For "eh" a character string specifies the base model for the respective statistical tests.
                Since "eh" searches for models between a defined minimal and a maximal model, a string vector of length 2 is another option.
                In this case a minimal and a maximal model is to be given.
                The model formula has to start with the first lowercase letters of the alphabet, e.g. "abc,cde". Variable names cannot be given.
                For "combined" the model cannot be given, a start model will be specified by \code{\link{gm.screening}}.
                 }
  \item{eh.strategy}{ See the strategy option of \code{eh} from \pkg{CoCo}.
                     }
  \item{criterion}{ Defines the strategy for edge selection: Likelihood-Ratio (default), AIC or BIC. Selection may be abbreviated.
                   }
  \item{tests}{ If TRUE, all accepted and rejected models are tested against each other. }
  \item{\dots}{ Additional options for the selection procedure. See \code{backward} from \pkg{CoCo} (forward) or \code{eh} from \pkg{CoCo}. }
}
\value{
A list containing:
  \item{"accepted" }{ String vector of the accepted model(s). }
  \item{"rejected" }{ String vector of the rejected model(s). }
  \item{"base"}{ Character string of model the selection strategy starts with. }
  \item{"strategy"}{ Character string indicating what type of selection strategy was performed. }
  \item{"tests"}{ Matrix containing the test results of the accepted vs. the rejected models 
                    (available if at least one selected model is submodel of another one).
                    Otherwise the accepted and rejected models are tested against the saturated model. }
  \item{"variable names"}{ Matrix that assigns a letter to each variable used in the model formula. }
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
\note{ CoCo may crash R occasionally on a bad day. We tried to prevent all crash causes for discrete data. 
        The R version of CoCo can only handle one CoCoObject.
        You may run "endCoCo()" if you use CoCo also outside this package or the function aborts.\cr
        \code{Criterion} equalling "bic" will work from CoCo version 0.1.7.5 on. \cr
        
        The EH-strategy is very time consuming. This depends on the number of variables, since
        the total number of possible models doubles with each additional variable.
}
\seealso{ \code{\link{gm.mim}}, \code{backward}, \code{eh}  from \pkg{CoCo}, \code{\link{gm.screening}} }
\examples{
\dontrun{
## until CoCo is available again ...
  data(wam)
  gm.coco(wam)
  ### giving many options to the strategy
  gm.coco(wam,criterion="b",recursive=TRUE,follow=TRUE,decomposable.mode=TRUE,
            coherent=TRUE)
  
  ### giving base
  gm.coco(wam,strategy="e",model=c("ab,cd","ae,be"))
  gm.coco(wam,strategy="f",model="abc,cd,de,f")
}
}
\keyword{ models }
\keyword{ multivariate }
