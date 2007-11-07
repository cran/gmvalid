\name{gm.boot.coco}
\alias{gm.boot.coco}
\title{ Graphical model validation using the bootstrap (CoCo). }
\description{
  Validates a discrete undirected graphical model using the bootstrap.
  To achieve that the CoCo package is used.
  Relative frequencies of the bootstrapped models, cliques or edges are counted.
}
\usage{
gm.boot.coco(N, data, strategy = c("backwards", "forwards", "eh"),
            calculations = c("subgraph", "diff", "edge", "clique"),
            model = FALSE, ...)
}
\arguments{
  \item{N}{ Number of bootstrap replications. }
  \item{data}{ Data frame or a table (array). Variables should have names, \code{data} has to be discrete. }
  \item{strategy}{ Type of model selection. "backwards" searches for not significant edges to delete, starting from the saturated model as default (\code{\link[CoCoRaw]{backward}}).
                    "forwards" adds significant edges, starting from the main effects model (\code{\link[CoCoRaw]{forward}}). The "Edwards-Havranek"
                    model search rejects complete models in every step and finishes with one or more accepted models (\code{\link[CoCoRaw]{eh}}).
                    The default strategy is "backwards". Selections may be abbreviated.
                     }
  \item{calculations}{ String vector specifying the analysis methods. "subgraph", "clique" and "edge" show the frequency of occurence 
                        in all selected models of the bootstrap samples, while "diff" counts the differences
                        in the bootstrap replications compared to the edges selected from the original data set.
                        The frequency of the selected models of all bootstrap samples are always calculated. By default all calculations are done.
                        Selections may be abbreviated.
                         }
  \item{model}{ Character string specifying a start model for "backwards" and "forwards" selection procedure.
                For "eh" a character string would specify the base model for tests.
                Since "eh" searches for models between a defined minimal and a maximal model, a string vector of length 2 is another option.
                In this case a minimal and a maximal model is to be given.
                The model formulae has to start with the first lowercase letters of the alphabet, e.g. "abc,cde". Variable names cannot be given.
                 }
  \item{\dots}{ Additional options to be given to the selection procedures. See \code{\link[CoCoRaw]{backward}} (forward) or \code{\link[CoCoRaw]{eh}}. }
}
\details{
  This function uses a nonparametric bootstrap.\cr
  For your information about the advancements of the bootstrap, some run numbers will be displayed.
}
\value{
A list containing:
  \item{"bootstrapped models" }{ Relative frequency vector of selected models with model formulas as categories. }
  \item{"bootstrapped subgraphs" }{ Relative frequencies vector of subgraphs that equals cliques or parts of a clique.
                                    Returned if calculation "subgraph" is selected.
                                     }
   \item{"bootstrapped cliques" }{ Relative frequency vector of selected cliques. Returned if calculation "clique" is selected.
                                     }
   \item{"bootstrapped edges"}{ Matrix with relative frequencies of the edges. Returned if calculation "edge" is selected. }
   \item{"original model"}{ Character string giving the selected model using the original unsampled data. Returned if calculation "diff" is selected. }
   \item{"differences from original data set"}{ List with frequencies of more, less and different edges.
                                                Sorted by occurence.
                                                Returned if calculation "diff" is selected. }
   \item{"variable names"}{ Matrix that assigns a letter to each variable that is used in the model formulas. }
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
\note{ CoCo may crash R occasionally on a bad day. The R version of CoCo can only handle one CoCoObject.
        You may run "endCoCo()" if you use CoCo also outside this package or the function aborts.\cr 
        Running this function may take a while. A bootstrap with
        data of size 10.000 times 8 and 1.000 replications may take the computer about 2 minutes.
}
\seealso{ \code{\link{gm.boot.mim}} }
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
}
\keyword{ nonparametric }
\keyword{ multivariate }
\keyword{ models }
