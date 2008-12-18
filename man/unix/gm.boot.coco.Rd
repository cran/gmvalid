\name{gm.boot.coco}
\alias{gm.boot.coco}
\title{ Graphical model validation using the bootstrap (CoCo). }
\description{
  Validates a discrete undirected graphical model using the bootstrap.
  Relative frequencies of the bootstrapped models, cliques or edges are counted.
  To achieve that the CoCo package is used.
}
\usage{
gm.boot.coco(N, data, strategy = c("backwards", "forwards", "combined"),
            calculations = c("diff", "edge", "clique"),
            model = FALSE, criterion = c("lr","aic","bic"), ...)
}
\arguments{
  \item{N}{ Number of bootstrap replications. }
  \item{data}{ Data frame or array. Variables need to be discrete and should have names. }
  \item{strategy}{ Type of model selection. "backwards" eliminates not significant edges, starting from the saturated model as default 
                    (\code{\link[CoCoRaw]{backward}}).
                    "forwards" adds significant edges, starting from the main effects model (\code{\link[CoCoRaw]{forward}}). 
                    The "combined" strategy is a 3 step procedure: \code{\link{gm.screening}}, "backwards" and then "forwards".
                    The default strategy is "backwards". Selections may be abbreviated.
                     }
  \item{calculations}{ String vector specifying the analysis methods; "clique" and "edge" show the frequency of occurence 
                        in all selected models of the bootstrap samples, while "diff" counts the edge differences
                        of the bootstrap replications compared to the edges selected from the original data set.
                        The frequency of the selected models of all bootstrap samples are always calculated. By default all calculations are done.
                        Selections may be abbreviated.
                         }
  \item{model}{ Character string specifying a start model for "backwards" or "forwards".
                The model formula has to start with the first lowercase letters of the alphabet, e.g. "abc,cde". Variable names cannot be given.
                For "combined" the model cannot be given, a start model will be specified by \code{\link{gm.screening}}.
                 }
  \item{criterion}{ Defines the strategy for edge selection: Likelihood-Ratio (default), AIC or BIC. Selection may be abbreviated.
                   }
  \item{\dots}{ Additional options for the selection procedure. See \code{\link[CoCoRaw]{backward}} (forward). }
}
\details{
  This function uses a nonparametric bootstrap.\cr
  For your information about the advancements of the bootstrap, some run numbers are displayed.
}
\value{
A list containing:
  \item{"bootstrapped models" }{ Matrix with the selected models in first column and their selection frequency in the second. }
   \item{"bootstrapped cliques" }{ Relative frequenciy vector of selected cliques. Returned if calculation "clique" is selected.
                                     }
   \item{"edge frequencies"}{ Matrix of relative edge frequencies. Returned if calculation "edge" is selected. }
   \item{"original model"}{ Character string of selected model using the original unsampled data. Returned if calculation "diff" is selected. }
   \item{"edge differences"}{ List of frequencies of more, less and by absolute value different edges (See argument \code{calculations}.
                                                Sorted by occurence.
                                                Returned if calculation "diff" is selected. }
  \item{"replications"}{ Number of bootstrap replications. }
  \item{"variable names"}{ Matrix that assigns a letter to each variable that is used in the model formulae. }
}
\references{ Efron B, Tibshirani RJ (1993) 
            \emph{An Introduction to the Bootstrap.}
             Chapman & Hall
             
               Badsberg JH (2001) 
               \emph{A guide to CoCo.}
  Journal of Statistical Software, 6(4).
             }
\author{ 
  Ronja Foraita, Fabian Sobotka \cr
  Bremen Institute for Prevention Research and Social Medicine \cr
  (BIPS)  \url{http://www.bips.uni-bremen.de}
     }
\note{ CoCo may crash R occasionally. The R version of CoCo can only handle one CoCoObject.
        You may run "endCoCo()" if you use CoCo also outside this package or the function aborts.\cr 
        \code{Criterion} equalling "bic" will work from CoCo version 0.1.7.5 on.
}
\seealso{ \code{\link{gm.screening}} }
\examples{
  ### should provide good results because of simulated data
  gm <- gm.modelsim(2000,"ABC,CDE")
  gm.boot.coco(50,gm,recursive=TRUE)
  
  ### on real data sets a forward bootstrap seems to have better results
  data(wynder)
  gm.boot.coco(100,wynder,strategy="f",calculations=c("s","e"),decomposable.mode=TRUE)
  
  ### with a given model
  data(wam)
  gm.boot.coco(10,wam,model="ab,bcde,cdef")
  ### by AIC
  gm.boot.coco(300,wam,strategy="f",criterion="a",recursive=TRUE,follow=TRUE,all.significant=FALSE)
}
\keyword{ nonparametric }
\keyword{ multivariate }
\keyword{ models }
