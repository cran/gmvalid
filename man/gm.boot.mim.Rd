\name{gm.boot.mim}
\alias{gm.boot.mim}
\title{ Graphical model validation using the bootstrap (MIM) }
\description{
  Validates a discrete undirected graphical model using the bootstrap.
  Relative frequencies of the bootstrapped models, cliques or edges are counted.
  Make sure that MIM is running.
}
\usage{
gm.boot.mim(N, data, strategy = c("backwards", "forwards", "eh", "combined"),
            calculations = c("diff", "edge", "clique"),
            model = FALSE, options = "")
}
\arguments{
  \item{N}{ Number of bootstrap replications. }
  \item{data}{ Data frame or array. Variables need to be discrete and should have names. }
  \item{strategy}{ Type of model selection. "backwards" eliminates not significant edges, starting from the saturated model as default.
                    "forwards" adds significant edges, starting from the main effects model. The "eh" ("Edwards-Havranek")
                    model search rejects complete models in every step and finishes with one or more accepted models.
                    The "combined" strategy is a 3 step procedure: \code{\link{gm.screening}}, "backwards" and then "forwards".
                    The default strategy is "backwards". Selections may be abbreviated.
                     }
  \item{calculations}{ String vector specifying the analysis methods; "clique" and "edge" show the frequency of occurence 
                        in all selected models of the bootstrap samples, while "diff" counts the edge differences
                        of the bootstrap replications compared to the edges selected from the original data set.
                        The frequency of the selected models of all bootstrap samples are always calculated. By default all calculations are done.
                        Selections may be abbreviated.
                     }
  \item{model}{ Character string to specify a start model for "backwards" or "forwards" selection procedures.
                For the "eh" procedure a minimal and a maximal model has to be assigned in one string
                connected with " - " (see Example). For "combined" the model cannot be given, a start model 
                will be specified by \code{\link{gm.screening}}.
                The model formula has to start with the first lowercase letters of the alphabet, e.g. "abc,cde". 
                Variable names cannot be given.
                 }
  \item{options}{ Optional character string to specify further options for the search strategy. 
                  Possible options can be found in the MIM help searching for
                    "stepwise" (backwards, forwards) or "startsearch" (eh).  See \code{details}.
                    }
}
\details{
  This function uses a nonparametric bootstrap.\cr
  For your information about the advancements of the bootstrap, some run numbers are displayed.
  
  MIM options for stepwise procedures (backwards, forwards):\cr
  "A" - uses the AIC as selection criterion\cr
  "B" - uses the BIC as selection criterion\cr
  "J" - joggles between backward and forward\cr
  "N" - non coherent mode\cr
  "U" - unrestricted, allows for non-decomposable models;
  
  MIM options for the eh modelsearch:\cr
  (positive number) - maximum number of models fitted\cr
  white space\cr
  (letters) - "U" for upward search, "D" for downward, default is both.
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

  Edwards D (2000) 
  \emph{An Introduction to Graphical Modelling.}
  Second Edition, Springer Verlag. 

             }
\author{ 
  Ronja Foraita, Fabian Sobotka \cr
  Bremen Institute for Prevention Research and Social Medicine \cr
  (BIPS)  \url{http://www.bips.uni-bremen.de}
     }
\note{ The function requires the MIM program.
        Make sure that it is running before using the function.\cr 
        The package mimR will only work properly if your R working directory path 
        does not contain hyphens ("-"). For mimR requirements we refer to the package's 
        help page \code{\link[mimR]{mimR}} and its homepage 
        \url{http://gbi.agrsci.dk/~shd/public/mimR/index.html}. mimR requires the Rgraphviz package. 
        Therefore you need to add "Bioconductor" to your R repositories.\cr
        
        The EH-strategy is very time consuming. This depends on the number of variables, since
        the total number of possible models doubles with each additional variable.
}
\seealso{ \code{\link{gm.boot.coco}}, \code{\link{gm.screening}} }
\examples{
  ### Examples work!
\dontrun{
  ### should provide good results because of simulated data
  gm.a <- gm.modelsim(2000,"ABC,CDE")
  gm.boot.mim(50,gm.a)
  
  ### on real data sets a forward bootstrap seems to have better results
  data(wynder)
  gm.boot.mim(100,wynder,strategy="f",calculations=c("s","e"),options="u")
  
  ### with model given
  data(wam)
  gm.boot.mim(10,wam,model="a,bcde,cdef")
  
  ### EH-strategy
  gm.boot.mim(50,wam,strategy="eh",model="a,bc,de,f - abcde,bcdef")
  }
}
\keyword{ nonparametric }
\keyword{ multivariate }
\keyword{ models }
