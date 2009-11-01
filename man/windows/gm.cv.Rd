\name{gm.cv}
\alias{gm.cv}
\title{ Cross Validation for Graphical (Chain) Models }
\description{
 Predicts a binary outcome variable in a given graphical (chain) model
 using k-fold cross validation.  
}
\usage{
gm.cv(k, data, outcome=1, strategy = c("backwards", "forwards"),
        chain = FALSE, options="", conf.level = 0.95, show.output = FALSE)
}
\arguments{
  \item{k}{ Number of folds the data is being split into to estimate the success of prediction in the cross validation. }
  \item{data}{ Data frame or array. Variables need to be discrete and should have names. }
  \item{outcome}{ Variable index of the outcome variable. Default is 1.}
  \item{strategy}{ Type of model selection. "backwards" searches for not significant edges to delete, starting from the saturated model as default.
                    "forwards" adds significant edges, starting from the main effects model. 
                    The default strategy is "backwards". Selections may be abbreviated.
                    }
  \item{chain}{ Character string to specify the block structure of directed graphical models.
                The syntax is "vs1|vs2|vs3" where vs1,vs2,vs3 are sets of variables and the variables in vs1 are prior to those in vs2 etc. Otherwise the model will be undirected. 
                Only lowercase letters are allowed!}
  \item{options}{ Character string to specify further options for the search strategy. Possible options can be found in the MIM help when searching for
                    "stepwise" (backwards, forwards) or "startsearch" (eh).  See \code{details}.
  }
  \item{conf.level}{ Confidence level of the interval (default is 0.95). }
  \item{show.output}{ Logical. logical. If TRUE then the output of the selection procedure is printed out of MIM. }
}
\details{
  Outcome variable must be the first in \code{data}, categorized as 
  1="unaffected" and 2="affected".

  The procedure is as follows
  
  1. Data Processing:\cr
     The data set is splitted into k folds.
  
  2. Model Selection:\cr
     A model is selected using the observations from k - j-th folds, j=1,...,k.
     P-values of edges are stored. 
     
  3. Calculate Risk:\cr
     A ratio table of being affected rather than unaffected is calulated for the joint probability of all influences associated with the outcome variable.
     This is done for every fold.\cr    
     If the ratio is greater or equal to one, the risk is set to 2, otherwise it is set to 1. For each fold a risk table is generated.
  
  4. Prediction:\cr
     The risk table from step(3), corresponding to the k-j-th fold, is used to predict the observations in the j-th fold.
     The prediction "PRED" is compared to the real outcome "OUT" in each fold by calculating the success probability using the following formula:
     
     success.prob = 1/n ( 1 - SUM(1:n)|OUT - PRED|),\cr
     n = number of oberservations in the j-th fold \cr
     
  MIM options for stepwise procedures (backwards, forwards):\cr
  "A" - uses the AIC as selection criterion\cr
  "B" - uses the BIC as selection criterion\cr
  "J" - joggles between backward and forward\cr
  "N" - non coherent mode\cr
  "U" - unrestricted, allows for non-decomposable models; 
}
\value{
A list containing:
  \item{ pvalue }{ Matrix of p-values for each edge in each fold. 
                  NA's mark missing edges. }
  \item{ ratio }{ Ratio table (see details, step (3)) for the best prediction model. }
  \item{ risk }{ Risk table (see details, step (3)) for the best prediction model. }
  \item{ success }{ Matrix containing the best prediction models of each fold, the number of edges that point
                   to the outcome variable and the probability of successful prediction (see details, step (4)). 
                   The initial blockstructure in the prediction step is: "variable set of influences | outcome variable".   
                   A given \code{chain} is only used during the model selection step.                   
                   The success probability is calculated using those clique structures the outcome variable is involved in.
                   }
  \item{ statistics }{ Matrix with Sensitivity, Specificity, Precision and Cohen's Kappa measures for the best
                        prediction model.
                        }
  \item{ variable names }{ Matrix that assigns a letter to each variable that is used in the model formulas. }
  
Additionally the graph of the best prediction model is plotted in MIM.
  
}
\references{ 
  Foraita R (2008) 
  \emph{Outcome prediction in graphical (chain) models using cross validation.}
  Slides. Please contact \email{foraita@bips.uni-bremen.de}.
  
  Edwards D (2000) 
  \emph{An Introduction to Graphical Modelling.}
  Second Edition, Springer Verlag. 
  }                  
\author{ 
  Ronja Foraita, Fabian Sobotka \cr
  Bremen Institute for Prevention Research and Social Medicine \cr
  (BIPS)  \url{http://www.bips.uni-bremen.de}
   }
\note{ 
The function requires the MIM program.
Make sure that it is running before using the function.
mimR will only work properly if your temporary directory has a path
where every folder has a name containing only 8 letters or less. mimR needs the Rgraphviz package.
Therefore you will have to add "Bioconductor" to your R repositories.
}
\seealso{ \code{\link{gm.mim}} }
\examples{
  ### Examples work under Windows (but not under Unix)
\dontrun{
  ABC <- gm.modelsim(500,"ABC,CD")
  out <- gm.cv(5,data=ABC, strategy="b", options="b")
  out
  
  ### DAG using a stepwise selection that joggles between backward and forward
  ### The direction of the chain is: d -> b -> c -> a
  out.dag <- gm.cv(3,data=ABC,option="j",chain="d|b|c|a")  
  
  ### Chain graph using BIC as selection criteria and allowing for 
  ### non-decomposable models
  cg <- gm.modelsim(1000,"ABD,BCE")  
  out.cg <- gm.cv(3,data=cg,option="bu",chain="cb|de|a")
  }
  ### Examples work neither under Windows nor Unix
  \dontrun{gm.cv(3,data=ABC,chain="DBD|A")  # you have to use lowercase letters
            gm.cv(3,data=ABC,chain="dca|b") # a is supposed to be outcome variable 
                                  # and thus have to be in the very right block    
            }
}
\keyword{ models }
\keyword{ regression }
