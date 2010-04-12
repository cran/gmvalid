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
to analyze discrete graphical models utilizing the CoCo package
and to quantify associations or interactions.
    
    Furthermore, the uncertainty of a selected graphical model can be described using the bootstrap or
the best prediction model can be evaluated for a dichotomous outcome variable and several discrete influences
using cross validation.
}
\details{
\tabular{ll}{
Package: \tab gmvalid\cr
Type: \tab Package\cr
Version: \tab 1.22\cr
Date: \tab 2010-04-12\cr
License: \tab GPL (>=2) \cr
}

\itemize{
\item This package requires the packages grid, epitools, gtools and CoCo.
\item If CoCo calculations fail, please type "endCoCo()" to remove temporary files. 
\item CoCo crashes occasionally. If you have recurring problems with a specific data set, we suggest to reboot the computer.
\item Model formulae have to be specified according to CoCo model formulae.
}

}
\author{
  Ronja Foraita, Fabian Sobotka \cr
  Bremen Institute for Prevention Research and Social Medicine \cr
  (BIPS)  \url{http://www.bips.uni-bremen.de}
}
\encoding{latin1}
\references{     
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
 \code{ \link[CoCo:CoCo-package]{CoCo} }
}
\examples{
    ### Generates a data frame given a dependence model
    gm.a <- gm.modelsim(1000,"ABC,CDE")
    
    ### Modelselection with graphical output
    gm.analysis(gm.a)   
    
    ### Model validation using the bootstrap 
    gm.boot.coco(100,gm.a,recursive=TRUE,follow=TRUE)

    ### Testing interaction on the penetrance scale
    ### using the conditional synergy index (CSI)
    gm.csi(1,2,3,data=gm.a)

    ### Testing interaction on a additivity scale
    ### using the synergy index (S)
    gm.si(1,2,3,data=gm.a)   

    ### Gamma Coefficient B indpendent D given C
    gm.gamma(2,4,data=gm.a,conditions=3)

}
