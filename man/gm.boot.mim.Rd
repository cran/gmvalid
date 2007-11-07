\name{gm.boot.mim}
\alias{gm.boot.mim}
\title{ Graphical model validation using the bootstrap (MIM). }
\description{
  Validates a discrete undirected graphical model using the bootstrap.
  To achieve that the program MIM is used.
  Relative frequencies of the bootstrapped models, cliques or edges are counted.
}
\usage{
gm.boot.mim(N, data, strategy = c("backwards", "forwards", "eh"),
            calculations = c("subgraph", "diff", "edge", "clique"),
            model = FALSE, options = "")
}
\arguments{
  \item{N}{ Number of bootstrap replications. }
  \item{data}{ Data frame or a table (array). Variables should have names, \code{data} has to be discrete. }
  \item{strategy}{ Type of model selection. "backwards" searches for not significant edges to delete, starting from the saturated model as default.
                    "forwards" adds significant edges, starting from the main effects model. The "Edwards-Havranek"
                    model search rejects complete models in every step and finishes with one or more accepted models.
                    The default strategy is "backwards". Selections may be abbreviated.
                     }
  \item{calculations}{ String vector specifying the analysis methods. "subgraph", "clique" and "edge" show the frequency of occurence 
                        in all selected models of the bootstrap samples, while "diff" counts the differences
                        in the bootstrap replications compared to the edges selected from the original data set.
                        The frequency of the selected models of all bootstrap samples are always calculated. By default all calculations are done.
                        Selections may be abbreviated.
                         }
  \item{model}{ Character string specifying a start model for "backwards" and "forwards" selection procedure.
                For the "eh" procedure a minimal and a maximal model has to be assigned in one string
                connected with a " - " (see Example). The model formula has to start with the first lowercase letters of the alphabet, e.g. "abc,cde". Variable names cannot be given.
                 }
  \item{options}{ Optional character string specifying further options for the search strategy. Possible options can be found in the MIM help searching for
                    "stepwise" (backwards, forwards) or "startsearch" (eh).  See \code{details}.
                    }
}
\details{
  This function uses a nonparametric bootstrap.\cr
  For your information about the advancements of the bootstrap, some run numbers will be displayed.
  
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
  \item{"bootstrapped models" }{ Relative frequency vector of selected models with model formulas as categories. }
  \item{"bootstrapped subgraphs" }{ Relative frequencies vector of subgraphs that equals cliques or parts of a clique.
                                    Returned if calculation "subgraph" is selected.
                                     }
   \item{"bootstrapped cliques" }{ Relative frequencies vector selected cliques. Returned if calculation "clique" is selected.
                                     }
   \item{"bootstrapped edges"}{ Matrix with relative frequencies of the edges. Returned if calculation "edge" is selected. }
   \item{"original model"}{ Character string giving the selected model using the original unsampled data. Returned if calculation "diff" is selected. }
   \item{"differences from original data set"}{ List with frequencies of more, less and different edges (See argument \code{calculations}.
                                                Sorted by occurence.
                                                Returned if calculation "diff" is selected. }
  \item{"variable names"}{ Matrix that assigns a letter to each variable that is used in the model formulas. }
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
        The package mimR will only work properly if your temporary directory has a path
        where every folder has a name containing only 8 letters or less. mimR needs the Rgraphviz package. Therefore you will have to add "Bioconductor" to your R repositories.\cr
        Running this function may take a while. A bootstrap with
        data of size 10.000 times 8 and 1.000 replications may take the computer about 4 minutes.
}
\seealso{ \code{\link{gm.boot.coco}} }
\examples{
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
\keyword{ nonparametric }
\keyword{ multivariate }
\keyword{ models }
