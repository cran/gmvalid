\name{gmvalid-package}
\alias{gmvalid-package}
\alias{gmvalid}
\docType{package}
\title{
    Validation of graphical models 
}
\description{
    This package provides functions among others that can be used to analyse
graphical models. This includes e.g. the possibility
to simulate data sets given a dependence model,
to analyze discrete graphical models utilizing the MIM program or the CoCo package
and to quantify associations or interactions.
    
    Furthermore, the uncertainty of a selected graphical model can be described using the bootstrap or
the best prediction model can be evaluated for a dichotomous outcome variable and several discrete influences
using cross validation.
}
\details{
\tabular{ll}{
Package: \tab gmvalid\cr
Type: \tab Package\cr
Version: \tab 1.23\cr
Date: \tab 2012-10-14\cr
License: \tab GPL (>=2) \cr
}

\itemize{
\item This package requires the packages grid, epitools, gtools and CoCo.
\item The package mimR is optional as well as the corresponding software MIM that can be downloaded at \url{http://www.hypergraph.dk}.
\item The mimR package itself requires the package Rgraphviz that is found in the Bioconductor repository.
\item Further the packages grid, epitools and gtools are required.
\item Before using any functions that use MIM make sure that the program is runnning.
\item The package mimR will only work properly if your R working directory path does not contain hyphens ("-").
\item For mimR requirements we refer to the package's help page \code{\link[mimR]{mimR}} and its homepage 
      \url{http://gbi.agrsci.dk/~shd/public/mimR/index.html}.
\item If CoCo calculations fail, please type "endCoCo()" to remove temporary files. 
\item CoCo crashes occasionally. If you have recurring problems with a specific data set, we suggest to use the MIM-option or to reboot the computer.
\item Model formulae have to be specified according to MIM model formulas or respectively CoCo model formulae.
\item To provide a package version that also runs on unix and therefore has no MIM functionalities, we had to put all MIM related examples in
      "Not Run". But for Windows the examples work.
}

}
\author{
  Ronja Foraita \cr
  Bremen Institute for Prevention Research and Social Medicine \cr
  (BIPS)  \url{http://www.bips.uni-bremen.de}
  
  Fabian Sobotka \cr
  Georg August University Goettingen \cr
  \url{http://www.uni-goettingen.de/de/350901.html}
}
\encoding{latin1}
\references{   
  > MIM (\url{http://www.hypergraph.dk/}) \cr
  Edwards D (2002) 
  \emph{An Introduction to Graphical Modelling.}
  Springer
  
  > mimR (\url{http://genetics.agrsci.dk/~sorenh/mimR/index.html})\cr
  Højsgaard S (2004) 
  \emph{The mimR package for graphical modelling in R.}
  Journal of Statistical Software, 11(6).
  
  > CoCo (\url{http://www.badsberg.eu})\cr
  Badsberg JH (2001) 
  \emph{A guide to CoCo.}
  Journal of Statistical Software, 6(4).
 
  > CSI \cr
  Foraita R (2008) 
  \emph{A conditional synergy index to assess biological interaction.}
  \url{http://nbn-resolving.de/urn:nbn:de:gbv:46-diss000111139}

}
\note{
  This work has been supported by the German Research Foundation\cr
  (DFG: \url{http://www.dfg.de})
  under grant scheme PI 345/2-1.
}
\keyword{ datagen }
\keyword{ models }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ graphs }
\keyword{ htest }
\seealso{
 \code{\link[mimR]{mimR}} and CoCo in the archive \url{http://cran.r-project.org/src/contrib/Archive/CoCo/}
}
\examples{
    ### Generates a data frame given a dependence model
    gm.a <- gm.modelsim(1000,"ABC,CDE")
\dontrun{
## until CoCo is available again ...
    ### Modelselection with graphical output
    gm.analysis(gm.a)   
    
    ### Model validation using the bootstrap 
    gm.boot.coco(100,gm.a,recursive=TRUE,follow=TRUE)

    ### Model prediction using cross validation
    gm.cv(3,data=gm.a,strategy="f",options="b")
}
    ### Testing interaction on the penetrance scale
    ### using the conditional synergy index (CSI)
    gm.csi(1,2,3,data=gm.a)

    ### Testing interaction on a additivity scale
    ### using the synergy index (S)
    gm.si(1,2,3,data=gm.a)   

    ### Gamma Coefficient B indpendent D given C
    gm.gamma(2,4,data=gm.a,conditions=3)

}
