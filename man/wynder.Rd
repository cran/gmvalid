\name{wynder}
\alias{wynder}
\docType{data}
\title{ Alcohol, Smoking and Oral Cancer }
\description{
  Alcohol drinking and smoking among cases of oral cancer and controls.
}
\usage{data(wynder)}
\format{
  A data frame with 443 observations on the following 3 variables.
  \describe{
    \item{\code{Alc}}{a numeric vector indicating the units of alcohol.}
    \item{\code{Smoking}}{a numeric vector indicating the units of smoked cigarettes.}
    \item{\code{Group}}{a numeric vector indicating whether the subject is "control" (1) or "case" (2).}
  }
}
\details{
  One unit of alcohol equals 1 oz (30 cc) of whiskey or 8 oz (237 cc) of beer.
}
\source{
  Wynder EE, Bross IJ, Feldman R (1957) \emph{A study of etiological factors in cancer of the mouth.}
  Cancer, 10:1300-1323.
}
\references{
  Rothman K (1974) \emph{The estimation of synergy or antagonism.}
  American Journal of Epidemiology, 103(5):506-511
}
\examples{
  data(wynder)
  tab <- table(wynder)
  dimnames(tab) <- list(c("<1","1-6",">6"), c("<15","16-34",">34"), c("controls","cases"))  
  names(dimnames(tab)) <- c("Alcohol (unit/day)","Smoking (cigarettes/day)","Group")  
  tab
}
\keyword{datasets}
