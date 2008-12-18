\name{gm.mim}
\alias{gm.mim}
\title{ Graphical model selection (MIM) }
\description{
  Selects one or more graphical models (depending on strategy) using the program MIM.
  Different types of model selection are available. A base model can be optionally defined.
  }
\usage{
gm.mim(data, strategy = c("backwards", "forwards", "eh", "combined"),
        model = FALSE, chain = FALSE, options = "", tests = FALSE)
}
\arguments{
  \item{data}{ Data frame or array. Variables need to be discrete and should have names. }
  \item{strategy}{ Type of model selection. "backwards" eliminates not significant edges, starting from the saturated model as default.
                    "forwards" adds significant edges, starting from the main effects model. The "eh" ("Edwards-Havranek")
                    model search rejects complete models in every step and finishes with one or more accepted models.
                    The "combined" strategy is a 3 step procedure: \code{\link{gm.screening}}, "backwards" and then "forwards".
                    The default strategy is "backwards". Selections may be abbreviated.
                     }
  \item{model}{ Character string to specify a start model for "backwards" or "forwards" selection procedures.
                For the "eh" procedure a minimal and a maximal model has to be assigned in one string
                connected with " - " (see Example). For "combined" the model cannot be given, a start model 
                will be specified by \code{\link{gm.screening}}.
                The model formula has to start with the first lowercase letters of the alphabet, e.g. "abc,cde". 
                Variable names cannot be given.
                 }
  \item{chain}{ Character string to specify the block structure of directed graphical models.
                The syntax is "vs1|vs2|vs3" where vs1,vs2,vs3 are sets of variables and the variables in vs1 are prior to those in vs2 etc.
                The formula has to be given with the first lowercase letters of the alphabet, e.g. "a|bc|d". Variable names cannot be given.
                For more information search the MIM help for "block overview". Does not work with "combined" selection strategy.
                 }
  \item{options}{  Character string specifying further options for the search strategy. Possible options can be found in the MIM help when searching for
                    "stepwise" (backwards, forwards) or "startsearch" (eh).  See \code{details}.
                    }
  \item{tests}{ If TRUE, all accepted and rejected models are tested against the saturated model, p-values to all edges are returned. }
}
\details{
  MIM options for stepwise procedures (backwards, forwards):\cr
  "A" - uses the AIC as selection criterion\cr
  "B" - uses the BIC as selection criterion\cr
  "J" - joggles between backward and forward\cr
  "N" - non coherent mode\cr
  "U" - unrestricted, allows for non-decomposable models;
  
  MIM options for the eh modelsearch:\cr
  (positive number) - maximum number of models fitted\cr
  white space\cr
  (letters) - "U" for upward search, "D" for downward, default is alternating.
  
}
\value{
A list containing:
  \item{"accepted" }{ String vector of accepted model(s). }
  \item{"rejected" }{ String vector of rejected model(s). }
  \item{"base"}{ Character string of model the selection starts with. NA for "eh". }
  \item{"strategy"}{ Character string indicating what type of selection strategy was performed. }
  \item{"tests against saturated"}{ Matrix giving the test results of the accepted and the rejected models against the saturated model.
                                    Only returned if \code{tests} is TRUE.
   }
  \item{"p values"}{ Upper-tri matrix of p-values of those edges in accepted models. NA if edge is misseing.
                    Only returned if \code{tests} is TRUE.
                     }
  \item{"variable names"}{ Matrix that assigns a letter to each variable that is used in the model formulas. }
}
\references{ 
  Edwards D (2000) 
  \emph{An Introduction to Graphical Modelling.}
  Second Edition, Springer Verlag. 
   }
\author{ 
  Ronja Foraita, Fabian Sobotka \cr
  Bremen Institute for Prevention Research and Social Medicine \cr
  (BIPS)  \url{http://www.bips.uni-bremen.de}
   }
\note{  The function requires the MIM program.
        Make sure that it is running before using the function. The package mimR will only work properly if your R working directory path 
        does not contain hyphens ("-"). For mimR requirements we refer to the package's help page \code{\link[mimR:mimR-package]{mimR}} and its homepage 
        \url{http://gbi.agrsci.dk/~shd/public/mimR/index.html}. mimR needs the Rgraphviz package.
        Therefore you will have to add "Bioconductor" to your R repositories.
        
        The EH-strategy is very time consuming. This depends on the number of variables, since
        the total number of possible models doubles with each additional variable.
}
\seealso{ \code{\link{gm.coco}} }
\examples{
  ### Examples work!
\dontrun{
  data(wam)
  gm.mim(wam)
  ### giving strategy
  gm.mim(wam,strategy="e")
  
  ### giving minimal and maximal model
  gm.mim(wam,strategy="e",model="a,bc,de,f - abcd,cdef")
  ### giving block structure
  gm.mim(wam,strategy="f",model="a,b,c|abc,de|abcd,ef",chain="abc|de|f",options="BNU")
  }
}
\keyword{ models }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
